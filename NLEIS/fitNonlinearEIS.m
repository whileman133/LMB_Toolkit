function data = fitNonlinearEIS(labSpectra,labLinearFit,varargin)
%FITNONLINEAREIS Regress nonlinear EIS model to lab impedance spectra.
%
% -- Usage --
% regressionData = fitNonlinearEIS(labSpectra,labLinearFit) regresses
%   the reduced-layer Warburg-resistance second-harmonic impedance model 
%   for LMB to the second-harmonoc impedance spectra in the structure
%   labSpectra created with loadLabNLEIS. labLinearFit is the data
%   structure from the linear EIS regression.
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
%   .arg         Structure of arguments supplied to the function.
%
% -- Changelog --
% 2023.07.17 | Created | Wesley Hileman <whileman@uccs.edu>

parser = inputParser;
parser.addRequired('labSpectra',@isstruct);
parser.addRequired('labLinearFit',@isstruct);
parser.addParameter('WeightFcn',[],@(x)isa(x,'function_handle'));
parser.parse(labSpectra,labLinearFit,varargin{:});
arg = parser.Results;  % structure of validated arguments

% Constants.
TdegC = arg.labSpectra.TdegC;
T = TdegC+273.15;

% Fetch OCP parameters fit to laboratory data.
ocpmodel = MSMR(labLinearFit.values.pos);

% Fetch linear impedance measured in the laboratory.
freqLab = labSpectra.h2.freq;
Zlab = labSpectra.h2.Z;

% Compute true SOC and lithiation for each SOC setpoint.
zmin = ocpmodel.zmin;
zmax = ocpmodel.zmax;
socPctTrue = labLinearFit.socPctTrue;
thetaTrue = zmax + (socPctTrue/100)*(zmin-zmax);

% Precompute OCP parameters at each SOC setpoint.
ocpData = ocpmodel.ocp('theta',thetaTrue,'TdegC',TdegC);

% Build model -------------------------------------------------------------
v = labLinearFit.values;

% Known parameters.
% OCP.
params.const.Q = fastopt.param('fix',v.const.Q);
params.pos.theta0 = fastopt.param('fix',v.const.theta0);
params.pos.theta100 = fastopt.param('fix',v.const.theta100);
params.pos.X = fastopt.param('fix',v.pos.X);
params.pos.U0 = fastopt.param('fix',v.pos.U0);
params.pos.omega = fastopt.param('fix',v.pos.omega);
% Electrolyte.
params.const.psi = fastopt.param('fix',v.const.psi);
params.const.W = fastopt.param('fix',v.const.W);
params.pos.tauW = fastopt.param('fix',v.pos.tauW);
params.pos.kappa = fastopt.param('fix',v.pos.kappa);
params.eff.tauW = fastopt.param('fix',v.pos.tauW);
params.eff.kappa = fastopt.param('fix',v.pos.kappa);
% Porous electrode.
params.pos.Rdl = fastopt.param('fix',v.pos.Rdl);
params.pos.Cdl = fastopt.param('fix',v.pos.Cdl);
params.pos.nDL = fastopt.param('fix',v.pos.nDL);
params.pos.Rf = fastopt.param('fix',v.pos.Rf);
params.pos.k0SplineTheta = fastopt.param('fix',v.pos.k0SplineTheta);
params.pos.k0Spline = fastopt.param('fix',v.pos.k0Spline);
params.pos.DsSplineTheta = fastopt.param('fix',v.pos.DsSplineTheta);
params.pos.DsSpline = fastopt.param('fix',v.pos.DsSpline);
params.pos.nF = fastopt.param('fix',v.pos.nF);
params.pos.sigma = fastopt.param('fix',v.pos.sigma);
% Lithium-metal electrode.
params.neg.Rdl = fastopt.param('fix',v.neg.Rdl);
params.neg.Cdl = fastopt.param('fix',v.neg.Cdl);
params.neg.nDL = fastopt.param('fix',v.neg.nDL);
params.neg.Rf = fastopt.param('fix',v.neg.Rf);
params.neg.k0 = fastopt.param('fix',v.neg.k0);

% Free paramters.
% Reaction symmetry factors.
params.neg.alpha = fastopt.param;
params.pos.alphaSpline = fastopt.param('len',ocpmodel.J);

modelspec = fastopt.modelspec(params);


% Define optimization bounds ----------------------------------------------
lb.neg.alpha = 0;                         ub.pos.alpha = 1;
lb.pos.alphaSpline = zeros(ocpmodel.J,1); ub.pos.alphaSpline = ones(ocpmodel.J,1);

init.neg.alpha = 0.5;
init.pos.alphaSpline = 0.5*ones(ocpmodel.J,1);

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
data.arg = arg;
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