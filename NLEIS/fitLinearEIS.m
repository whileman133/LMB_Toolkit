function data = fitLinearEIS(labSpectra,labOCPFit,initialModel,varargin)
%FITLINEAREIS Regress linear EIS model to laboratory impedance spectra.
%
% -- Usage --
% regressionData = fitLinearEIS(labSpectra,labOCPFit,initialModel) regresses
%   the reduced-layer Warburg-resistance transfer-function model for LMB 
%   to the linear impedance spectra in the structure labSpectra created 
%   with loadLabNLEIS. labOCPFit is the OCP model to use in performing the 
%   EIS regression. initialModel is a model containing initial estimates of 
%   parameter values for the cell.
%
% Output structure:
%   .estimate    Cell model containing estimated parameter values
%   .modelspec   fastopt specification for EIS model
%   .values      Structure of estimated parameter values
%   .initial     Structure of initial parameter values
%   .lb          Structure of lower bounds for parameters
%   .ub          Structure of upper bounds for parameters
%   .Zmodel      Matrix of model-predicted impedance: d1=frequency d2=soc
%   .Zlab        Matrix of lab impedance: d1=frequency d2=soc
%   .freq        Cyclic frequency vector
%   .socPctTrue  SOC vector, calculated from QdisAh
%   .TdegC       Temperature
%
% -- Changelog --
% 2023.06.30 | Created | Wesley Hileman <whileman@uccs.edu>

parser = inputParser;
parser.addRequired('labSpectra',@isstruct);
parser.addRequired('labOCPFit',@isstruct);
parser.addRequired('initialModel',@isstruct);
parser.addParameter('WeightFcn',[],@(x)isa(x,'function_handle'));
parser.parse(labSpectra,labOCPFit,initialModel,varargin{:});
arg = parser.Results;  % structure of validated arguments

% Constants.
tplus0 = 0.4;    % Li+ transference number, guess for estimating psi
R = TB.const.R;  % molar gas constant [J/mol/K]
F = TB.const.F;  % Faraday's constant [C/mol]
T = arg.labSpectra.TdegC+273.15; % absolute temperature [K]

% Fetch OCP parameters fit to laboratory data.
ocpmodel = MSMR(labOCPFit.MSMR);
ocptest = labOCPFit.ocptest;

% Fetch linear impedance measured in the laboratory.
freqLab = labSpectra.lin.freq;
Zlab = labSpectra.lin.Z;
TdegC = labSpectra.TdegC;

% Compute true SOC and lithiation for each SOC setpoint.
zmin = ocpmodel.zmin;
zmax = ocpmodel.zmax;
QdisAhCum = cumsum(labSpectra.QdisAh);
socPctTrue = 100*(1-QdisAhCum/ocptest.QAh);
thetaTrue = zmax + (socPctTrue/100)*(zmin-zmax);

% Precompute OCP parameters at each SOC setpoint.
ocpData = ocpmodel.ocp('theta',thetaTrue,'TdegC',TdegC);

% Load OCP estimate into initial model.
% Convert initial model to reduced-layer Warburg-resistance model.
% Fetch initial values of model parameters.
p = struct;
p.pos.U0 = ocpmodel.Uj0;
p.pos.X = ocpmodel.Xj;
p.pos.omega = ocpmodel.Wj;
p.pos.theta0 = ocpmodel.zmax;
p.pos.theta100 = ocpmodel.zmin;
initialModel = setCellParam(initialModel,p);
initialModel = convertCellModel(initialModel,'RLSWRM');
initial = getCellParams(initialModel,'TdegC',labSpectra.TdegC);

% Build model -------------------------------------------------------------
% Known parameters.
params.const.Q = fastopt.param('fix',ocptest.QAh);
params.pos.theta0 = fastopt.param('fix',ocpmodel.zmax);
params.pos.theta100 = fastopt.param('fix',ocpmodel.zmin);
params.pos.X = fastopt.param('fix',ocpmodel.Xj);
params.pos.U0 = fastopt.param('fix',ocpmodel.Uj0);
params.pos.omega = fastopt.param('fix',ocpmodel.Wj);

% Fixed parameters.
% Rdl(p|n) are low for the Sion cells, so fix Rdl(p|n)=0.
params.pos.Rdl = fastopt.param('fix',0);
params.neg.Rdl = fastopt.param('fix',0);
params.neg.Rf = fastopt.param('fix',0);  % Lump into tab resistance.
% Symmetry factors for the positive electrode are not very identifiable,
% we'll asume alpha=0.5 for all galleries.
params.pos.alpha = fastopt.param('fix',0.5*ones(ocpmodel.J,1));
% Symmetry factor for the negative electrode is not identifiable due to
% linear nature of the small-signal impedance model.
params.neg.alpha = fastopt.param('fix',0.5);
% Solid conductance does not influence impedance significatly.
params.pos.sigma = fastopt.param('fix',initial.pos.sigma);
% Fix lithiation of spline points.
params.pos.k0SplineTheta = fastopt.param('fix',initial.pos.k0SplineTheta);
params.pos.DsSplineTheta = fastopt.param('fix',initial.pos.DsSplineTheta);

% Electrolyte parameters.
% psi,W,tauW are not separately identifiable, so we fix psi=R/(F*(1-t+0)), 
% the value predicted by the Einstien relationship.
params.const.psi = fastopt.param('fix',R/F/(1-tplus0));
params.const.W = fastopt.param;
params.pos.tauW = fastopt.param('logscale',true);
params.pos.kappa = fastopt.param('logscale',true);
params.eff.tauW = fastopt.param('logscale',true);
params.eff.kappa = fastopt.param('logscale',true);

% Porous electrode parameters.
params.pos.DsSpline = fastopt.param('len',ocpmodel.J,'logscale',true);
params.pos.nF = fastopt.param;
params.pos.k0Spline = fastopt.param('len',ocpmodel.J,'logscale',true);
params.pos.Cdl = fastopt.param;
params.pos.nDL = fastopt.param;
params.pos.Rf = fastopt.param;

% Lithium-metal electrode parameters.
params.neg.k0 = fastopt.param('logscale',true);
params.neg.Cdl = fastopt.param;
params.neg.nDL = fastopt.param;

% Cell package parameters.
params.pkg.R0 = fastopt.param;  % Tab resistance.
params.pkg.L0 = fastopt.param;  % Package/cable inductace.

modelspec = fastopt.modelspec(params);


% Define optimization bounds ----------------------------------------------
% Initial guess; pack/unpack to set values of fixed parameters.
init = fastopt.unpack(fastopt.pack(initial,modelspec),modelspec);

% Electrolyte parameters.
lb.const.W = 0.1;                   ub.const.W = 10;
lb.eff.tauW = init.eff.tauW/100;    ub.eff.tauW = init.eff.tauW*100;
lb.pos.tauW = init.pos.tauW/100;    ub.pos.tauW = init.pos.tauW*100;
lb.eff.kappa = init.eff.kappa/100;  ub.eff.kappa = init.eff.kappa*100;
lb.pos.kappa = init.pos.kappa/100;  ub.pos.kappa = init.pos.kappa*100;

% Porous-electrode parameters.
lb.pos.DsSpline = 1e-9*ones(ocpmodel.J,1);ub.pos.DsSpline = 1*ones(ocpmodel.J,1);
lb.pos.nF = 0.1;                    ub.pos.nF = 1;
lb.pos.k0Spline = 1e-8*ones(ocpmodel.J,1);ub.pos.k0Spline = 1e6*ones(ocpmodel.J,1);
lb.pos.sigma = init.pos.sigma/10;   ub.pos.sigma = init.pos.sigma*10;
lb.pos.Cdl = init.pos.Cdl/10;       ub.pos.Cdl = init.pos.Cdl*100;
lb.pos.nDL = 0.5;                   ub.pos.nDL = 1;
lb.pos.Rf = init.pos.Rf/100;        ub.pos.Rf = init.pos.Rf*100;

% Lithium-metal electrode parameters.
lb.neg.k0 = init.neg.k0/1000;       ub.neg.k0 = init.neg.k0*1000;
lb.neg.Cdl = init.neg.Cdl/10;       ub.neg.Cdl = init.neg.Cdl*10;
lb.neg.nDL = 0.5;                   ub.neg.nDL = 1;

% Cell package parameters.
lb.pkg.R0 = 0;                      ub.pkg.R0 = init.pkg.R0*10;
lb.pkg.L0 = 0;                      ub.pkg.L0 = init.pkg.L0*100;

% Perform regression ------------------------------------------------------

% Collect weight matrix.
weights = ones(length(freqLab),length(socPctTrue));
if ~isempty(arg.WeightFcn)
    for idxSOC = 1:length(socPctTrue)
        for idxFreq = 1:length(freqLab)
            weights(idxFreq,idxSOC) = arg.WeightFcn( ...
                freqLab(idxFreq),socPctTrue(idxSOC));
        end % for
    end % for
end

[plotInit, plotUpdate] = uiImpedanceCallbacks( ...
    modelspec,Zlab,freqLab,socPctTrue,TdegC);
data = fastopt.uiparticleswarm(@cost,modelspec,init,lb,ub, ...
    'PlotInitializeFcn',plotInit,'PlotUpdateFcn',plotUpdate);

% Collect output data.
data.model = setCellParam(initialModel,data.values);
data.Zmodel = getLinearImpedance(data.values,freqLab,socPctTrue,TdegC,ocpData);
data.Zlab = Zlab;
data.freq = freqLab;
data.socPctTrue = socPctTrue;
data.TdegC = TdegC;
data.arg = struct;
data.arg.labSpectra = labSpectra;
data.arg.labOCPFit = labOCPFit;
data.arg.initialModel = initialModel;
data.type__ = 'ParameterEstimate';
data.origin__ = 'fitLinearEIS';

function J = cost(model)
    % Ignore singular matrix warnings in solving TF model (usu. occurs at high
    % frequencies). TODO: find high-frequency TF solution.
    warning('off','MATLAB:nearlySingularMatrix');

    % Calculate impedance predicted by the linear EIS model.
    Zmodel = getLinearImpedance(model,freqLab,socPctTrue,TdegC,ocpData);

    % Re-enable singular matrix warning.
    warning('on','MATLAB:nearlySingularMatrix');

    % Compute total residual between model impedance and measured
    % impedance across all spectra.
    J = sum(weights.*(abs(Zmodel-Zlab)./abs(Zlab)).^2,'all');
end % cost()

end % fitLinearEIS()