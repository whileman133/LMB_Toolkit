function fitData = fitRPT(meas,modelspec,init,lb,ub,varargin)
%FITRPT Update parameter values by regressing RPT model to measurements.
%  Particle-swarm method.
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
%   LB and UB are structures of lower and upper bounds for the
%     parameter values. 
%   INIT is a structure of starting parameter values.
%
% -- Data Format --
% Input: Measurement Structure (MEAS)
%  .halfCycle   : Structure containing processed half-cycle discharge data.
%    .socAvgPct : average cell SOC vector [%]
%    .vcell     : cell voltage vector [V]
%    .TdegC     : average cell temperature [degC]
%  .eis         : Structure containing processed (NL)EIS data.
%    .lin.Z     : matrix of linear spectra (dim1=freq, dim2=SOC) [V/A]
%    .lin.freq  : cyclic frequency vector for linear spectra [Hz]
%    .socPct    : SOC vector for both linear and second-harmonoc spectra [%]
%    .TdegC     : average cell temperature [degC]
%
% Output: Fit Data Structure (FITDATA)
%  .estimate        : structure of estimated parameter values
%  .predict.hcVcell : half-cycle voltage predicted by regressed model [V]
%  .predict.linZ    : linear spectra predicted by regressed model [V/A]
%  .arg             : structure of arguments supplied to the function
%
% -- Changelog --
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
    hcyc.ocpData = electrode.ocp('theta',hcyc.theta,'TdegC',hcyc.TdegC);
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

function J = cost(model)
    % Ignore singular matrix warnings in solving TF model (usu. occurs at high
    % frequencies). TODO: find high-frequency TF solution.
    warning('off','MATLAB:nearlySingularMatrix');

    % Calculate impedance predicted by the linear EIS model.
    Zmodel = getLinearImpedance( ...
        model,eis.,socPctTrue{km},TdegC(km),ocpData{km});
    % Compute total residual between model impedance and measured
    % impedance across all spectra.
    J = J + sum(weights{km}.*(abs(Zmodel-Zlab{km})./abs(Zlab{km})).^2,'all');

    % Re-enable singular matrix warning.
    warning('on','MATLAB:nearlySingularMatrix');
    
end

end