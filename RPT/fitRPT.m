function fitData = fitRPT(meas,modelspec,init,lb,ub,varargin)
%FITRPT Update parameter values by regressing RPT model to measurements
%  using particle-swarm optimization (PSO).
%
%  The Reference Performance Test (RPT) consists of a half-cycle discharge
%  followed by (NL)EIS at several SOC setpoints.
%
%  This function supports regression to both simulated measurements 
%  (for verification) and laboratory data (for application).
%
% -- Usage --
% fitData = fitRPT(meas,modelspec,init,lb,ub) regresses the RPT model to the
%   measurements in the structure MEAS using particle swarm optimization.
%   MODELSPEC is the fastopt.modelspec structure specifying the parameters 
%     over which to optimize the RPT model. 
%   LB and UB are structures of lower and upper bounds for the parameter 
%     values. 
%   INIT is a structure of starting parameter values.
%
% -- Data Format --
% Input: Measurement Structure (MEAS)
%  .halfCycle   : Structure containing processed half-cycle discharge data.
%    .socAvgPct : average cell SOC vector [%]
%    .Rdc       : dc resistance vector [Omega]
%    .TdegC     : average cell temperature [degC]
%  .eis         : Structure containing processed (NL)EIS data.
%    .lin.Z     : matrix of linear spectra (dim1=freq, dim2=SOC) [V/A]
%    .lin.freq  : cyclic frequency vector for linear spectra [Hz]
%    .socPct    : SOC vector [%]
%    .TdegC     : average cell temperature [degC]
%
% Output: Fit Data Structure (FITDATA)
%  .estimate        : structure of estimated parameter values
%  .lb              : structure of lower bounds (user could have changed)
%  .ub              : structure of upper bounds (user could have changed)
%  .J               : cost function value for the final estimate
%  .predict.Rdc     : dc resistance predicted by regressed model [V]
%  .predict.linZ    : linear spectra predicted by regressed model [V/A]
%  .arg             : structure of arguments supplied to this function
%
% -- Changelog --
% 2023.11.14 | Updated for new RPT protocol | Wes H.
% 2023.09.12 | Created | Wesley Hileman <whileman@uccs.edu>

parser = inputParser;
parser.addRequired('meas',@(x)isstruct(x)&&isscalar(x));
parser.addRequired('modelspec',@(x)isstruct(x)&&isscalar(x));
parser.addRequired('init',@(x)isstruct(x)&&isscalar(x));
parser.addRequired('lb',@(x)isstruct(x)&&isscalar(x));
parser.addRequired('ub',@(x)isstruct(x)&&isscalar(x));
parser.addParameter('EISWeightFcn',[],@(x)isa(x,'function_handle'));
parser.addParameter('UseParallel',true,@(x)islogical(x)&&isscalar(x));
parser.parse(meas,modelspec,init,lb,ub,varargin{:});
arg = parser.Results;  % struct of validated arguments

% Structures containing data for the respective tests.
hcyc = arg.meas.halfCycle;
eis = arg.meas.eis;

% Code optimization: precompute OCP when the OCP parameters are fixed.
fixed = fastopt.getFixedParams(arg.modelspec);
if all(isfield(fixed.pos,{'theta0','theta100','U0','X','omega'}))
    p = fixed.pos;
    electrode = MSMR(p);
    eis.theta = p.theta0+(eis.socPct/100)*(p.theta100-p.theta0);
    eis.ocpData = electrode.ocp('theta',eis.theta,'TdegC',eis.TdegC);
    hcyc.thetaAvg = p.theta0+(hcyc.socAvgPct/100)*(p.theta100-p.theta0);
    hcyc.ocpData = electrode.ocp('theta',hcyc.thetaAvg,'TdegC',hcyc.TdegC);
else
    eis.theta = [];
    eis.ocpData = [];
    hcyc.thetaAvg = [];
    hcyc.ocpData = [];
end

% Collect EIS weight matrix.
eis.weights = ones(length(eis.lin.freq),length(eis.socPct));
if ~isempty(arg.EISWeightFcn)
    for idxSOC = 1:length(eis.socPct)
        for idxFreq = 1:length(eis.lin.freq)
            eis.weights(idxFreq,idxSOC) = arg.EISWeightFcn( ...
                eis.lin.freq(idxFreq),eis.socPct(idxSOC),eis.TdegC);
        end % for
    end % for
end % if

% Perform regression.
[plotInit, plotUpdate] = uiRPTCallbacks(modelspec,hcyc,eis);
psoData = fastopt.uiparticleswarm(@cost,modelspec,init,lb,ub, ...
    'PlotInitializeFcn',plotInit,'PlotUpdateFcn',plotUpdate, ...
    'UseParallel',arg.UseParallel);

% Collect output.
fitData.estimate = psoData.values;
fitData.lb = psoData.lb;
fitData.ub = psoData.ub;
[fitData.J, fitData.predict.Rdc, fitData.predict.linZ] = ... 
    cost(psoData.values);
fitData.arg = arg;


function [J, Rmodel, Zmodel] = cost(model)
    J = 0;

    % 1. EIS component
    % Calculate impedance predicted by the linear EIS model.
    % Compute total residual between model impedance and measured
    % impedance across all spectra.
    warning('off','MATLAB:nearlySingularMatrix');
    Zmodel = getLinearImpedance( ...
        model,eis.lin.freq,eis.socPct,eis.TdegC,eis.ocpData);
    warning('on','MATLAB:nearlySingularMatrix');
    J = J + sum(eis.weights.*(abs(Zmodel-eis.lin.Z)./abs(eis.lin.Z)).^2,'all');
    
    % 2. Half-cycle component
    % Calculate Rdc predicted by perturbation model.
    % Compute total residual between model and dc resistance measured in 
    % the laboratory.
    if isempty(hcyc.thetaAvg)
        thetaAvg = model.pos.theta0+...
            (hcyc.socAvgPct/100)*(model.pos.theta100-model.pos.theta0);
    else
        thetaAvg = hcyc.thetaAvg;
    end
    dcData = getDcResistance( ...
        model,thetaAvg,'TdegC',hcyc.TdegC,'ocpData',hcyc.ocpData);
    Rmodel = dcData.Rdc;
    J = J + sum((Rmodel-hcyc.Rdc).^2./hcyc.Rdc.^2,'all');
end

end