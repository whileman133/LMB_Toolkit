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
% 2023.11.04 | Incorporate double-layer and CPE effects | Wes H.
% 2023.06.02 | Coerce output into column vector | Wesley H.
% 2022.08.22 | Created | Wesley Hileman <whileman@uccs.edu>

parser = inputParser;
parser.addRequired('model',@isstruct);
parser.addRequired('thetaAvg',@(x)isnumeric(x)&&isvector(x));
parser.addParameter('TdegC',25,@(x)isnumeric(x)&&isscalar(x));
parser.addParameter('ocpData',[],@istruct);
parser.addParameter('ComputeRctj',false,@islogical)
parser.parse(cellModel,theta,varargin{:});
arg = parser.Results; % structure of validated arguments

% Ensure lithiation is a row vector.
theta = theta(:)';

% Setup physical constants.
f = TB.const.f(arg.TdegC); % (F/T) at TdegC.

% Fetch struct of cell parameters.
cellParams = getCellParams(cellModel,'TdegC',arg.TdegC);
% Fetch OCP, charge-transfer resistance, Ds at each lithiation point.
if ~isempty(arg.ocpData)
    % First choice: Use Uocv provided to function (cached vector, fastest).
    ocpData = arg.ocpData;
else
    % Second resort: compute the OCP using the MSMR parameters 
    % (will call fzero twice and interp1 once, even slower).
    msmr = MSMR(cellParams.pos);
    ocpData = msmr.ocp('theta',theta,'TdegC',arg.TdegC);
end
ctData = msmr.RctCachedOCP(cellParams.pos,ocpData);
dsData = msmr.DsCachedOCP(cellParams.pos,ocpData);

% Compute parameters needed to evalulate dc resistance.
Q = cellParams.const.Q;
W = cellParams.const.W;
Rc = cellParams.const.Rc;  % tab resistance
theta0p = cellParams.pos.theta0;
theta100p = cellParams.pos.theta100;
nFp = cellParams.pos.nF;
tauFp = cellParams.pos.tauF;
nDLp = cellParams.pos.nDL;
tauDLp = cellParams.pos.tauDL;
sp = cellParams.pos.sigma;
kp = cellParams.pos.kappa;
Rfp = cellParams.pos.Rf;
Rdlp = cellParams.pos.Rdl;
Cdlp = cellParams.pos.Cdl;
dUocpp = ocpData.dUocp; % pos electrode OCP slope (vector)
Rctp = ctData.Rct;      % pos electrode charge-transfer resistance (vector)
Dsp = dsData.Ds;        % pos electrode diffusion coefficient (vector)
Dsdcp = Dsp.^nFp*tauFp.^(nFp-1); % diffusion coefficient @ dc (vector)
Cdldcp = Cdlp.^nDLp.*tauDLp.^(1-nDLp); % double-layer cap @ dc
Csp = -3600*Q./dUocpp./abs(theta100p-theta0p); % cell cap @ dc
Rfn = cellParams.neg.Rf;
k0n = cellParams.neg.k0;
Rctn = 1./f./k0n;
if isfield(cellParams,'eff')
    % eff layer combines dll and sep layers
    kd = Inf;
    ks = cellParams.eff.kappa;
else
    % individual dll and sep layers
    kd = cellParams.dll.kappa;
    ks = cellParams.sep.kappa;
end

% Compute SOC-independent resistance.
R0 = Rc + Rctn + Rfn + (1+W)*(1/kd+1/ks+1/3/kp) + 1/3/sp;

% Compute diffusion resistance.
Rd = 1./15./Csp./Dsdcp;

% Finally, calculate the total dc resistance.
Rdc = R0+( ... 
    (Rdlp+Rfp).*Cdldcp.^2 + 2.*Rfp.*Cdldcp.*Csp + (Rctp+Rd+Rfp).*Csp.^2 ...
    - tauDLp.*(1-nDLp).*Cdldcp - tauFp.*(1-nFp).*Csp ...
)./(Cdldcp+Csp).^2;

% Assign output struct.
data.Rdc = Rdc;
data.Uocp = ctData.Uocp;
data.parts.R0 = R0(:);
data.parts.Rd = Rd(:);
data.parts.Rctp = Rctp(:);
data.parts.Rctn = Rctn(:);
data.arg = arg;
data.origin__ = 'getDcResistance';

end
