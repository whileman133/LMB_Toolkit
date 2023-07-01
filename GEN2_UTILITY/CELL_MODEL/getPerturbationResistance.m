function data = getPerturbationResistance(model, thetaAvg, varargin)
%GETPERTURBATIONRESISTANCE Compute Baker-Verbrugge perturbation resistance
%  of a LMB cell. Lumped-parameter model.
%
% Rtotal = GETPERTURBATIONRESISTANCE(model,thetaAvg) computes the
%   perturbation resistance of the LMB cell specified by MODEL at the average 
%   lithiation values given in the vector THETAAVG. RTOTAL is a column
%   vector of the perturbation resistances.
%
% [Rtotal,parts] = GETPERTURBATIONRESISTANCE(...) also returns
%   a structure PARTS of individual resistances making up the perturbation 
%   resistance:
%      parts.R0 = SOC-invariant part (equivalent series resistance) [scalar].
%      parts.Rct_n = charge-transfer resistance of the negative electrode,
%         a component of R0 [scalar].
%      parts.Rdiff = resistance due to SOC-dependent solid diffusion [vector].
%      parts.Rct_p = total charge-transfer resistance of the intercalate
%        electrode interface, SOC variant [vector].
%   If the optional argument pair ('computeRctj',true) is also supplied,
%   then the PARTS structure has an additional element:
%      parts.Rctj_p = charge-transfer resistance associated with each MSMR
%        gallery [matrix]. Columns correspond to SOC setpoints. The total
%        charge-transfer resistance is given by: Rct_p = sum(Rctj_p).
%
% [Rtotal,parts,U] =  GETPERTURBATIONRESISTANCE(...) also returns the OCV 
%   evalulated at each lithiation setpoint THETAAVG (i.e., the dynamic 
%   equilibrium solution). U is computed internally in the course of 
%   determining the  MSMR gallery partial lithiations xj, which are needed 
%   to compute Rctp.
%
% [...] = GETPERTURBATIONRESISTANCE(...,'TdegC',T) performs the calculation
%   at temperature T instead of the default 25degC.
%
% -- Performance options --
% [...] = GETPERTURBATIONRESISTANCE(...,'Uocv',Uocv) performs the
%   calculation using the OCV vector UOCV instead of computing the OCV
%   for each lithiation point in THETAAVG. This can speed up the computation 
%   if this function is called repeatedly inside an optimization routine
%   each time with the same THETAAVG, temperature, and MSMR OCP parameters.
%   If the MSMR parameters change each iteration, this option cannot be used.
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
% [1] Daniel R. Baker and Mark W. Verbrugge 2021 J. Electrochem. Soc. 168 050526
%
% -- Changelog --
% 2023.06.02 | Coerce output into column vector | Wesley Hileman
% 2022.08.22 | Created | Wesley Hileman <whileman@uccs.edu>

parser = inputParser;
parser.addRequired('model',@isstruct);
parser.addRequired('thetaAvg',@(x)isnumeric(x)&&isvector(x));
parser.addParameter('TdegC',25,@(x)isnumeric(x)&&isscalar(x));
parser.addParameter('Uocv',[],@isnumeric);
parser.addParameter('ComputeRctj',false,@islogical)
parser.parse(model,thetaAvg,varargin{:});
p = parser.Results; % structure of validated parameter values

T = p.TdegC+273.15;
f = TB.const.F/TB.const.R/T;
Uocv = p.Uocv;
computeRctj = p.ComputeRctj;

% Ensure lithiation is a row vector.
thetaAvg = thetaAvg(:)';

if isfield(model,'function')
    % Toolbox cell model.
    isReg = @(reg)isfield(model.function,reg);
    getParam = @(reg,p)model.function.(reg).(p)(0,T);
else
    % Set of model parameters already evalulated at setpoint.
    isReg = @(reg)isfield(model,reg);
    getParam = @(reg,p)model.(reg).(p);
end

% Compute series resistance component (does not vary with SOC).
if isfield(model.const,'R0')
    % Model specifies R0 directly.
    Rct_n = NaN;
    R0 = model.const.R0;
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
        kappa_eff = getParam('eff','kappa');
        R0 = (Rf_p+Rct_n+Rf_n) + 1/sigma_p/3 + ...
             (1+W)*(1/kappa_p/3 + 1/kappa_eff);
    else
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
Rdiff = abs(theta100-theta0)/f/Q/Dsref./thetaAvg./(1-thetaAvg)/5/10800;

% Calculate Rct(pos) at each stoichiometry setpoint.
U0 = getParam('pos','U0');          U0 = U0(:);
X = getParam('pos','X');            X = X(:);
omega = getParam('pos','omega');    omega = omega(:);
alpha = getParam('pos','alpha');    alpha = alpha(:);
k0_p = getParam('pos','k0');        k0_p = k0_p(:);
if ~isempty(Uocv)
    % First choice: Use Uocv provided to function (cached vector, fastest).
    U = Uocv;
elseif isfield(model,'function')
    % Next best: Fetch OCP from toolbox model (calls interp1, slightly slower).
    U = model.function.pos.Uocp(thetaAvg,T);
else
    % Last resort: compute the OCP using the MSMR parameters 
    % (will call fzero twice and interp1 once, even slower).
    msmr = MSMR(model.pos);
    U = msmr.ocp('theta',thetaAvg,'TdegC',p.TdegC);
end
U = U(:)';
xj = X./(1+exp(f*(U-U0)./omega));
Rct_p = 1/f./sum(k0_p.*(xj).^(omega.*alpha).*(X-xj).^(omega.*(1-alpha)));
if computeRctj
    Rctj_p = 1/f./(k0_p.*(xj).^(omega.*alpha).*(X-xj).^(omega.*(1-alpha)));
end

% Finally, calculate the perturbation resistance.
Rtotal = R0 + Rct_p + Rdiff;
Rtotal = Rtotal(:);  % Force into column vector.

% Assign individual components of the resistance.
parts.R0 = R0(:);
parts.Rct_n = Rct_n(:);
parts.Rdiff = Rdiff(:);
parts.Rct_p = Rct_p(:);
if computeRctj
    parts.Rctj_p = Rctj_p;
end

data.Rtotal = Rtotal;
data.parts = parts;
data.U = U(:);
data.param = p;
data.origin__ = 'getPerturbationResistance';

end
