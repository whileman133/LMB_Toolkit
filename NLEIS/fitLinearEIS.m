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
% Supply datasets at multiple temperatures by specifying LABSPECTRA as a 
% structure array, containing a scalar structure for each temperature 
% dataset.
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
% 2023.09.19 | Fix error computing theta, wrong theta min/max | Wes H.
% 2023.08.27 | Store separate IS for Ds/i0 based on distance | Wes H.
% 2023.08.10 | Add option to save intermediate solutions (IS) | Wes H.
% 2023.07.26 | Update for multiple temperatures | Wes H.
% 2023.07.24 | Update for multiple diffusion and kinetics models | Wes H.
% 2023.06.30 | Created | Wesley Hileman <whileman@uccs.edu>

isdstype = @(x)any(strcmpi(x,{'msmr','linear','spline'}));
isi0type = @(x)any(strcmpi(x,{'linear','spline'}));
parser = inputParser;
parser.addRequired('labSpectra',@(x)isstruct(x)&&isvector(x));
parser.addRequired('labOCPFit',@(x)isstruct(x)&&isscalar(x));
parser.addRequired('initialModel',@isstruct);
parser.addParameter('WeightFcn',[],@(x)isa(x,'function_handle'));
parser.addParameter('SolidDiffusionModel','spline',isdstype);
parser.addParameter('KineticsModel','spline',isi0type);
parser.addParameter('TrefdegC',25,@(x)isnumeric(x)&&isscalar(x));
% NOTE: IntermediateSoln stuff is disabled when UseParallel option is true;
% MATLAB can't easily share data-structures between threads.
parser.addParameter('IntermediateSolns',100,@(x)isnumeric(x)&&isscalar(x)); % num previous minimizers to keep
parser.addParameter('IntermediateSolnEpsilon',struct('Ds',0.03,'k0',0.03),@(x)isscalar(x)&&istruct(x));
parser.addParameter('UseParallel',true,@(x)islogical(x)&&isscalar(x));
parser.parse(labSpectra,labOCPFit,initialModel,varargin{:});
arg = parser.Results;  % structure of validated arguments

% Determine the multiplicity of the dataset, i.e. the number of
% temperatures included.
multiplicity = length(arg.labSpectra);

% Constants.
tplus0 = 0.4;    % Li+ transference number, guess for estimating psi
R = TB.const.R;  % molar gas constant [J/mol/K]
F = TB.const.F;  % Faraday's constant [C/mol]

% Fetch OCP parameters fit to laboratory data.
ocpmodel = MSMR(labOCPFit.MSMR);
ocptest = labOCPFit.ocptest;
% !!! Important: We need to compute theta min/max from the voltage limits of 
% the cell, NOT the values in labOCPFit.MSMR, as our lab data may have been
% augmented with additional data collected over wider lithiation range.
ocpTmp = ocpmodel.ocp('voltage',[ocptest.vmin ocptest.vmax]);
thetamin = min(ocpTmp.theta);
thetamax = max(ocpTmp.theta);

% Collect linear impedance measured in the laboratory.
freqLab = cell(multiplicity,1);
Zlab = cell(multiplicity,1);
for m = 1:multiplicity
    freqLab{m} = arg.labSpectra(m).lin.freq;
    Zlab{m} = arg.labSpectra(m).lin.Z;
end
TdegC = [arg.labSpectra.TdegC].';  % temperature vector [degC].

% Compute true SOC and lithiation for each SOC setpoint.
% Precompute OCP parameters at each SOC setpoint.
socPctTrue = cell(multiplicity,1);
thetaTrue = cell(multiplicity,1);
ocpData = cell(multiplicity,1);
for m = 1:multiplicity
    QdisAhCum = cumsum(arg.labSpectra(m).QdisAh);
    socPctTrue{m} = 100*(1-QdisAhCum/ocptest.QAh);
    thetaTrue{m} = thetamax + (socPctTrue{m}/100)*(thetamin-thetamax);
    ocpData{m} = ocpmodel.ocp('theta',thetaTrue{m},'TdegC',TdegC(m));
end
nsocavg = round(mean(cellfun(@length,socPctTrue)));
minTheta = min(cellfun(@min,thetaTrue));
maxTheta = max(cellfun(@max,thetaTrue));

% Load OCP estimate into initial model.
% Convert initial model to reduced-layer Warburg-resistance model.
% Fetch initial values of model parameters.
p = struct;
p.pos.U0 = ocpmodel.Uj0;
p.pos.X = ocpmodel.Xj;
p.pos.omega = ocpmodel.Wj;
p.pos.theta0 = thetamax;
p.pos.theta100 = thetamin;
initialModel = setCellParam(initialModel,p);

% Convert cell model to reduced-layer Warburg-resistance model.
% Also convert the kinetics and diffusion sub-models.
thetaLin = linspace(minTheta,maxTheta,nsocavg);
if strcmpi(arg.KineticsModel,'linear')
    kStruct.type = 'linear';
    kStruct.theta = thetaLin(:);  % data point near each lab SOC setpoint
else
    kStruct.type = arg.KineticsModel;
end
if strcmpi(arg.SolidDiffusionModel,'linear')
    dStruct.type = 'linear';
    dStruct.theta = thetaLin(:);  % data point near each lab SOC setpoint
else
    dStruct.type = arg.SolidDiffusionModel;
end
initialModel = convertCellModel(initialModel,'RLWRM', ...
    'KineticsModel',kStruct,'SolidDiffusionModel',dStruct);
initial = getCellParams(initialModel,'TdegC',arg.TrefdegC);

% Build model -------------------------------------------------------------
% Known parameters.
params.const.Q = fastopt.param('fix',ocptest.QAh);
params.pos.theta0 = fastopt.param('fix',thetamax);
params.pos.theta100 = fastopt.param('fix',thetamin);
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
if strcmpi(arg.KineticsModel,'linear')
    params.pos.alphaLinear = fastopt.param('fix',0.5*ones(length(thetaLin),1));
elseif strcmpi(arg.KineticsModel,'spline')
    params.pos.alphaSpline = fastopt.param('fix',0.5*ones(ocpmodel.J,1));
else
    params.pos.alpha = fastopt.param('fix',0.5*ones(ocpmodel.J,1));
end
% Symmetry factor for the negative electrode is not identifiable due to
% linear nature of the small-signal impedance model.
params.neg.alpha = fastopt.param('fix',0.5);
% Solid conductance does not influence impedance significatly.
params.pos.sigma = fastopt.param('fix',initial.pos.sigma);
% Fix lithiation of interpolation points.
if any(strcmpi(arg.KineticsModel,{'linear','spline'}))
    params.pos.k0Theta = fastopt.param('fix',initial.pos.k0Theta);
end
if any(strcmpi(arg.SolidDiffusionModel,{'linear','spline'}))
    params.pos.DsTheta = fastopt.param('fix',initial.pos.DsTheta);
end
% psi,W,tauW are not separately identifiable, so we fix psi=R/(F*(1-t+0)), 
% the value predicted by the Einstien relationship.
params.const.psi = fastopt.param('fix',R/F/(1-tplus0));

% Electrolyte parameters.
params.const.W = fastopt.param('logscale',true);
params.pos.tauW = fastopt.param('logscale',true,'tempfcn','Eact','tempcoeff','-');
params.pos.kappa = fastopt.param('logscale',true,'tempfcn','Eact');
params.eff.tauW = fastopt.param('logscale',true,'tempfcn','Eact','tempcoeff','-');
params.eff.kappa = fastopt.param('logscale',true,'tempfcn','Eact');

% Porous electrode parameters.
if strcmpi(arg.SolidDiffusionModel,'msmr')
    params.pos.Dsref = fastopt.param('logscale',true,'tempfcn','Eact');
    params.pos.mD = fastopt.param;
elseif strcmpi(arg.SolidDiffusionModel,'linear')
    params.pos.DsLinear = fastopt.param('len',length(thetaLin),'logscale',true,'tempfcn','Eact');
else
    params.pos.DsSpline = fastopt.param('len',ocpmodel.J,'logscale',true,'tempfcn','Eact');
end
if strcmpi(arg.KineticsModel,'linear')
    params.pos.k0Linear = fastopt.param('len',length(thetaLin),'logscale',true,'tempfcn','Eact');
else
    params.pos.k0Spline = fastopt.param('len',ocpmodel.J,'logscale',true,'tempfcn','Eact');
end
params.pos.nF = fastopt.param;
params.pos.tauF = fastopt.param('logscale',true);
params.pos.Cdl = fastopt.param;
params.pos.nDL = fastopt.param;
params.pos.wDL = fastopt.param('logscale',true);
params.pos.Rf = fastopt.param;

% Lithium-metal electrode parameters.
params.neg.k0 = fastopt.param('logscale',true,'tempfcn','Eact');
params.neg.Cdl = fastopt.param;
params.neg.nDL = fastopt.param;
params.neg.wDL = fastopt.param('logscale',true);

% Cell package parameters.
params.pkg.R0 = fastopt.param('tempfcn','lut');  % Tab resistance.
params.pkg.L0 = fastopt.param('tempfcn','lut');  % Package/cable inductace.

modelspec = fastopt.modelspec(params, ...
    'tempsdegC',TdegC,'TrefdegC',arg.TrefdegC);


% Define optimization bounds ----------------------------------------------
% Initial guess; pack/unpack to set values of fixed parameters.
% Add activation energy parameters.
initial.pos.tauW_Eact = 1000;
initial.pos.kappa_Eact = 1000;
initial.eff.tauW_Eact = 1000;
initial.eff.kappa_Eact = 1000;
if strcmpi(arg.SolidDiffusionModel,'msmr')
    initial.pos.Dsref_Eact = 1000;
elseif strcmpi(arg.SolidDiffusionModel,'linear')
    initial.pos.DsLinear_Eact = 1000;
else
    initial.pos.DsSpline_Eact = 1000;
end
if strcmpi(arg.KineticsModel,'linear')
    initial.pos.k0Linear_Eact = 1000;
else
    initial.pos.k0Spline_Eact = 1000;
end
initial.neg.k0_Eact = 1000;
init = fastopt.unpack(fastopt.pack(initial,modelspec),modelspec);

% Electrolyte parameters.
lb.const.W = 0.2;                   ub.const.W = 200;
lb.eff.tauW = init.eff.tauW/100;    ub.eff.tauW = init.eff.tauW*100;
lb.eff.tauW_Eact = 0;               ub.eff.tauW_Eact = 100000;
lb.pos.tauW = init.pos.tauW/100;    ub.pos.tauW = init.pos.tauW*100;
lb.pos.tauW_Eact = 0;               ub.pos.tauW_Eact = 100000;
lb.eff.kappa = init.eff.kappa/100;  ub.eff.kappa = init.eff.kappa*100;
lb.eff.kappa_Eact = 0;              ub.eff.kappa_Eact = 100000;
lb.pos.kappa = init.pos.kappa/100;  ub.pos.kappa = init.pos.kappa*100;
lb.pos.kappa_Eact = 0;              ub.pos.kappa_Eact = 100000;

% Porous-electrode parameters.
if strcmpi(arg.SolidDiffusionModel,'msmr')
    lb.pos.Dsref = 1e-9;            ub.pos.Dsref = 1e3;
    lb.pos.Dsref_Eact = 0;          ub.pos.Dsref_Eact = 100000;
    lb.pos.mD = 1;                  ub.pos.mD = 5;
elseif strcmpi(arg.SolidDiffusionModel,'linear')
    lb.pos.DsLinear = 1e-5*ones(length(thetaLin),1);ub.pos.DsLinear = 1*ones(length(thetaLin),1);
    lb.pos.DsLinear_Eact = 0;       ub.pos.DsLinear_Eact = 100000;
else
    lb.pos.DsSpline = 1e-5*ones(ocpmodel.J,1);ub.pos.DsSpline = 1*ones(ocpmodel.J,1);
    lb.pos.DsSpline_Eact = 0;       ub.pos.DsSpline_Eact = 100000;
end
lb.pos.nF = 0.1;                    ub.pos.nF = 1;
lb.pos.tauF = 0.1;                  ub.pos.tauF = 10000;
if strcmpi(arg.KineticsModel,'linear')
    lb.pos.k0Linear = 1e-4*ones(length(thetaLin),1);ub.pos.k0Linear = 100*ones(length(thetaLin),1);
    lb.pos.k0Linear_Eact = 0;       ub.pos.k0Linear_Eact = 100000;
else
    lb.pos.k0Spline = 1e-4*ones(ocpmodel.J,1);ub.pos.k0Spline = 100*ones(ocpmodel.J,1);
    lb.pos.k0Spline_Eact = 0;       ub.pos.k0Spline_Eact = 100000;
end
lb.pos.sigma = init.pos.sigma/1000; ub.pos.sigma = init.pos.sigma*1000;
lb.pos.Cdl = init.pos.Cdl/10;       ub.pos.Cdl = init.pos.Cdl*100;
lb.pos.nDL = 0.5;                   ub.pos.nDL = 1;
lb.pos.wDL = 1e-6;                  ub.pos.wDL = 10;
lb.pos.Rf = init.pos.Rf/100;        ub.pos.Rf = init.pos.Rf*100;

% Lithium-metal electrode parameters.
lb.neg.k0 = init.neg.k0/1000;       ub.neg.k0 = init.neg.k0*1000;
lb.neg.k0_Eact = 0;                 ub.neg.k0_Eact = 100000;
lb.neg.Cdl = init.neg.Cdl/10;       ub.neg.Cdl = init.neg.Cdl*10;
lb.neg.nDL = 0.5;                   ub.neg.nDL = 1;
lb.neg.wDL = 1e-6;                  ub.neg.wDL = 10;

% Cell package parameters.
lb.pkg.R0 = 0;                      ub.pkg.R0 = init.pkg.R0*10;
lb.pkg.L0 = 0;                      ub.pkg.L0 = init.pkg.L0*100;


% Perform regression ------------------------------------------------------

% Collect weight matrix.
weights = cell(multiplicity,1);
for m = 1:multiplicity
    weights{m} = ones(length(freqLab{m}),length(socPctTrue{m}));
    if ~isempty(arg.WeightFcn)
        for idxSOC = 1:length(socPctTrue{m})
            for idxFreq = 1:length(freqLab{m})
                weights{m}(idxFreq,idxSOC) = arg.WeightFcn( ...
                    freqLab{m}(idxFreq),socPctTrue{m}(idxSOC),TdegC(m));
            end % for
        end % for
    end % if
end % for

% Allocate storage for intermediate solutions.
topJ.Ds = inf(1,arg.IntermediateSolns);
topJ.k0 = inf(1,arg.IntermediateSolns);
topSoln.Ds(arg.IntermediateSolns) = init;
topSoln.k0(arg.IntermediateSolns) = init;

% Perform regression.
[plotInit, plotUpdate] = uiImpedanceCallbacks( ...
    modelspec,Zlab,freqLab,socPctTrue,TdegC);
data = fastopt.uiparticleswarm(@cost,modelspec,init,lb,ub, ...
    'PlotInitializeFcn',plotInit,'PlotUpdateFcn',plotUpdate, ...
    'UseParallel',arg.UseParallel);

% Collect output data.
data.topJ = topJ;
data.topSoln = topSoln;
data.modelspec = modelspec;
models = fastopt.splittemps(data.values,modelspec);
for m = multiplicity:-1:1
    data.Zmodel{m} = getLinearImpedance( ...
        models(m),freqLab{m},socPctTrue{m},TdegC(m),ocpData{m});
end
data.Zlab = Zlab;
data.freq = freqLab;
data.socPctTrue = socPctTrue;
data.TdegC = TdegC;
data.initialModel = initialModel;
data.argpso = data.arg;
data = rmfield(data,{'origin__','arg'});
data.arg = arg;
data.type__ = 'ParameterEstimate';
data.origin__ = 'fitLinearEIS';

function J = cost(model)
    modelVect = fastopt.splittemps(model,modelspec);

    % Ignore singular matrix warnings in solving TF model (usu. occurs at high
    % frequencies). TODO: find high-frequency TF solution.
    warning('off','MATLAB:nearlySingularMatrix');

    J = 0;
    for km = 1:multiplicity
        % Calculate impedance predicted by the linear EIS model.
        Zmodel = getLinearImpedance( ...
            modelVect(km),freqLab{km},socPctTrue{km},TdegC(km),ocpData{km});
        % Compute total residual between model impedance and measured
        % impedance across all spectra.
        J = J + sum(weights{km}.*(abs(Zmodel-Zlab{km})./abs(Zlab{km})).^2,'all');
    end

    % Re-enable singular matrix warning.
    warning('on','MATLAB:nearlySingularMatrix');

    % Store intermediate solutions.
    if ~arg.UseParallel
        Jlt = J<topJ.Ds; 
        if any(Jlt) && strcmpi(arg.SolidDiffusionModel,'linear')
            if isinf(topJ.Ds(1))
                store = true;
            else
                pos = [topSoln.Ds(~isinf(topJ.Ds)).pos];
                n0Ds = lognormalize( ...
                    [pos.DsLinear],lb.pos.DsLinear,ub.pos.DsLinear);
                nDs = lognormalize( ...
                    model.pos.DsLinear,lb.pos.DsLinear,ub.pos.DsLinear);
                store = all(vecnorm(nDs-n0Ds)>arg.IntermediateSolnEpsilon.Ds);
            end % else
            if store
                ind = find(Jlt,1,'first');
                topJ.Ds = [topJ.Ds(1:ind-1) J topJ.Ds(ind:end-1)];
                topSoln.Ds = [topSoln.Ds(1:ind-1) model topSoln.Ds(ind:end-1)];
            end % if
        end % if
        Jlt = J<topJ.k0; 
        if any(Jlt) && strcmpi(arg.KineticsModel,'linear')
            if isinf(topJ.k0(1))
                store = true;
            else
                pos = [topSoln.k0(~isinf(topJ.k0)).pos];
                n0k0 = lognormalize( ...
                    [pos.k0Linear],lb.pos.k0Linear,ub.pos.k0Linear);
                nk0 = lognormalize( ...
                    model.pos.k0Linear,lb.pos.k0Linear,ub.pos.k0Linear);
                store = all(vecnorm(nk0-n0k0)>arg.IntermediateSolnEpsilon.k0);
            end % else
            if store
                ind = find(Jlt,1,'first');
                topJ.k0 = [topJ.k0(1:ind-1) J topJ.k0(ind:end-1)];
                topSoln.k0 = [topSoln.k0(1:ind-1) model topSoln.k0(ind:end-1)];
            end % if
        end % if
    end % if

    % Fragments from old cost function definitions:
    % 
    % 1. Phase-angle cost function
    % DeltaAngle = angle(Zmodel)-angle(Zlab);
    % indWrap = abs(DeltaAngle)>pi;
    % DeltaAngle(indWrap) = sign(DeltaAngle(indWrap)).*2*pi-abs(DeltaAngle(indWrap));
    % DeltaAngle = 180*DeltaAngle/pi;  % convert rad to degrees
    % NormAngle = 20;
    % J = sum(weights.*(DeltaAngle.^2./NormAngle.^2),'all');
end % cost()

function n = lognormalize(x,lb,ub)
    n = (log10(x)-log10(lb))./(log10(ub)-log10(lb));
end

end % fitLinearEIS()