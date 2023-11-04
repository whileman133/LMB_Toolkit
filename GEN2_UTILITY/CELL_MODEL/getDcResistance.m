function data = getDcResistance(cellModel, theta, varargin)
%GETDCRESISTANCE Compute low-frequency resistance of an LMB cell.
%
% data = GETDCRESISTANCE(cellModel,theta) computes the dc
%   perturbation resistance of the LMB cell specified by CELLMODEL at the 
%   lithiation values given in the vector THETA. The output, DATA, is a 
%   structure with the following fields:
%
%   .Rtotal       = column vector of the dc resistance [Ohm]
%   .Uocp         = column vector of cel OCP at THETA (i.e., 
%                   the dynamic equilibrium solution). U is computed 
%                   internally in the course of determining the MSMR 
%                   gallery partial lithiations xj, which are needed 
%                   to compute Rctp.
%   .parts.R0     = SOC-invariant part (equivalent series resistance) [scalar].
%   .parts.Rctn   = charge-transfer resistance of the negative electrode,
%                   a component of R0 [scalar].
%   .parts.Rd     = resistance due to SOC-dependent solid diffusion [vector].
%   .parts.Rctp   = total charge-transfer resistance of the intercalate
%                   electrode interface, SOC variant [vector].
%   .parts.Rctjp  = charge-transfer resistance associated with each MSMR
%                   gallery [matrix]. Columns correspond to SOC setpoints. 
%                   The total charge-transfer resistance is given by: 
%                   Rct_p = sum(Rctj_p).
%                   **Included only if the optional argument pair 
%                    ('ComputeRctj',true) is also supplied!**
%
% [...] = GETDCRESISTANCE(...,'TdegC',T) performs the calculation
%   at temperature T instead of the default 25degC.
%
% -- Performance options --
% [...] = GETDCRESISTANCE(...,'ocpData',ocpData) performs the
%   calculation using the OCV vector UOCV instead of computing the OCV
%   for each lithiation point in THETA. This can speed up the computation 
%   if this function is called repeatedly inside an optimization routine
%   each time with the same THETA, temperature, and MSMR OCP parameters.
%   If the MSMR parameters change each iteration, this option cannot be 
%   used.
%
% -- Background --
% The reduced-order perturbation approximation developed by 
% Baker and Verbrugge [1] models cell voltage as:
%   vcell(t) = Uocv(thetaAvg(t)) - iapp(t)*Rtotal(thetaAvg(t))
% where Rtotal(theta) is the SOC-dependent perturbation resistance of the
% cell. The approximation is valid provided iapp(t) appears constant on the
% scale of the solid diffusion time Rs^2/Ds. Uocv(thetaAvg(t)) is the
% dynamic-equilibrium solution. thetaAvg(t) is the average lithiation of
% the positive electrode ("absolute" SOC).
%
% This model also coincides with the dc limit of the transfer-function
% model developed at UCCS, hence the name "dc resistance".
%
% -- References --
% [1] Daniel R. Baker and Mark W. Verbrugge 2021 J. Electrochem. Soc. 168 050526
%
% -- Changelog --
% 2023.11.04 | 
% 2023.06.02 | Coerce output into column vector | Wesley Hileman
% 2022.08.22 | Created | Wesley Hileman <whileman@uccs.edu>

parser = inputParser;
parser.addRequired('model',@isstruct);
parser.addRequired('thetaAvg',@(x)isnumeric(x)&&isvector(x));
parser.addParameter('TdegC',25,@(x)isnumeric(x)&&isscalar(x));
parser.addParameter('ocpData',[],@istruct);
parser.addParameter('ComputeRctj',false,@islogical)
parser.parse(cellModel,theta,varargin{:});
arg = parser.Results; % structure of validated arguments

if isCellModel(cellModel)
    % Covert to legacy lumped-parameter model for use with code below.
    cellModel = convertCellModel(cellModel,'LLPM');
else
    % Assume a structure of parameter values was supplied instead;
    % no need to convert for code below.
end

T = arg.TdegC+273.15;
f = TB.const.F/TB.const.R/T;
computeRctj = arg.ComputeRctj;

% Ensure lithiation is a row vector.
theta = theta(:)';

% Define getters depending on the form of the model (functions or structure
% of values).
if isfield(cellModel,'function')
    % Toolbox cell model.
    isReg = @(reg)isfield(cellModel.function,reg);
    getReg = @(reg)cellModel.function.(reg);
    getParam = @(reg,p)cellModel.function.(reg).(p)(0,T);
else
    % Set of model parameters already evalulated at setpoint.
    isReg = @(reg)isfield(cellModel,reg);
    getReg = @(reg)cellModel.(reg);
    getParam = @(reg,p)cellModel.(reg).(p);
end

% Compute series resistance component (does not vary with SOC).
if isfield(cellModel.const,'R0')
    % Model specifies R0 directly.
    Rct_n = NaN;
    R0 = cellModel.const.R0;
else
    % Compute R0 from model parameters.
    W = getParam('const','W');
    Rf_p = getParam('pos','Rf');
    kappa_p = getParam('pos','kappa');
    sigma_p = getParam('pos','sigma');
    k0_n = getParam('neg','k0');
    Rf_n = getParam('neg','Rf');
    Rct_n = 1/f/k0_n;
    if isReg('eff')
        % eff layer combines dll and sep
        kappa_eff = getParam('eff','kappa');
        R0 = (Rf_p+Rct_n+Rf_n) + 1/sigma_p/3 + ...
             (1+W)*(1/kappa_p/3 + 1/kappa_eff);
    else
        % individual dll and sep layers
        kappa_s = getParam('sep','kappa');
        kappa_d = getParam('dll','kappa');
        R0 = (Rf_p+Rct_n+Rf_n) + 1/sigma_p/3 + ...
             (1+W)*(1/kappa_p/3 + 1/kappa_s + 1/kappa_d);
    end
end

% Calculate the diffusion resistance at each stoichiometry setpoint.
Q = getParam('const','Q');
theta0 = getParam('pos','theta0');
theta100 = getParam('pos','theta100');
Dsref = getParam('pos','Dsref');
Rd = abs(theta100-theta0)/f/Q/Dsref./theta./(1-theta)/5/10800;

% Calculate Rct(pos) at each stoichiometry setpoint.
if ~isempty(arg.ocpData)
    % First choice: Use Uocv provided to function (cached vector, fastest).
    ocpData = arg.ocpData;
else
    % Last resort: compute the OCP using the MSMR parameters 
    % (will call fzero twice and interp1 once, even slower).
    msmr = MSMR(getReg('pos'));
    ocpData = msmr.ocp('theta',theta,'TdegC',arg.TdegC);
end
ctData = msmr.RctCachedOCP(getReg('pos'),ocpData);
Rctp = ctData.Rct;
if computeRctj
    Rctjp = ctData.Rctj;
end

% Finally, calculate the perturbation resistance.
Rtotal = R0 + Rctp(:) + Rd(:);

% Assign individual components of the resistance.
parts.R0 = R0(:);
parts.Rctn = Rct_n(:);
parts.Rd = Rd(:);
parts.Rctp = Rctp(:);
if computeRctj
    parts.Rctjp = Rctjp;
end

data.Rtotal = Rtotal;
data.parts = parts;
data.U = ctData.Uocp;
data.arg = arg;
data.origin__ = 'getDcResistance';

end
