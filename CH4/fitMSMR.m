function data = fitMSMR(ocp, J, varargin)
%FITMSMR Regress an MSMR model to laboratory OCP measurements.
%
% * Uses ga and fmincon to perform the nonlinear optimization.
% * Ensures the galleries do not "swap" position so that lower and upper
%   bounds on U0,X,omega may be pinned to each gallery.
% * Ensures sum(Xj)=1.
%
% -- Usage --
% data = FITMSMR(ocp,J)
% data = FITMSMR(...,'initial',init,'eps',eps)
% data = FITMSMR(...,'lb',lb,'ub',lb)
% data = FITMSMR(...,'initial',init,'lb',lb,'ub',lb)
% data = FITMSMR(...,'fix',fix)
% data = FITMSMR(...,'Usep',Usep)
% data = FITMSMR(...,'w',w)
%
% -- Input --
% ocp   = struct or struct array of OCP data (see fields below)
% J     = number of galleries to fit
% vmin  = minimum voltage for fitting MSMR model [V] 
% vmax  = maximum voltage for fitting MSMR model [V] 
% init  = structure of initial values for regression (see fields below)
% eps   = use to automatically generate lower and upper bounds
%         from initial values: lb = (1-eps)*init, ub = (1+eps)*init
%         scalar or structure; scalar applies uniform boundaries to each 
%         parameter; structure for nonuniform boundaries (see fields below)
% lb    = struct of lower bounds (see fields below)
% ub    = struct of upper bounds (see fields below)
% fix   = structure of parameter values to fix to specific values and not
%         optimize over (see fields below)
% Usep  = minimum required separation of gallery potentials U0j [V]
% w     = relative weight placed on agreement of differential capacity 
%         vs OCP in the cost function [-]
%
% ocp is a struct or struct array with the following fields:
%
%   .TdegC     Temperature [degC]
%   .Uocp      OCP vector [V]
%   .Z         Relative composition vector [-]
%   .dZ        Differential capacity vector [1/V]
%
% init, eps, lb, ub, fix are structs with the following fields:
%
%   .U0        Column vector of gallery potentials U0j [V]
%   .X         Column vector of gallery partial compositions [-]
%   .omega     Column vector of omega_j values [-]
%   .thetamin  Minimim composition of electrode [-]
%   .thetamax  Maximum composition of electrode [-]
%
% -- Output --
% The output, DATA, is a struct with the following fields:
%
%   .est       the best-fit MSMR model found
%   .Jcost     cost associated with the best fit
%
% -- Changelog --
% 2024.01.21 | Cleaned up | Wesley Hileman <whileman@uccs.edu>

% Constants.
msmrFields.X.len = J;
msmrFields.U0.len = J;
msmrFields.omega.len = J;
msmrFields.thetamin.len = 1;
msmrFields.thetamax.len = 1;
msmrFieldnames = fieldnames(msmrFields);
msmrFieldnamesComposition = {'X','thetamin','thetamax'};
R = 8.3144598;      % Molar gas constant [J/mol K]
F = 96485.3329;     % Faraday constant [C/mol]

% Parse arguments.
parser = inputParser;
parser.addParameter('initial',[],@(x)isstruct(x)||isa(x,'MSMR')||isempty(x));
parser.addParameter('eps',0.5,@(x)isstruct(x)||isscalar(x)&&x>=0);
parser.addParameter('lb',struct,@isstruct);
parser.addParameter('ub',struct,@isstruct);
parser.addParameter('fix',struct,@isstruct);
parser.addParameter('Usep',0,@(x)x>=0);
parser.addParameter('w',0.001,@(x)x>=0);
parser.addParameter('gaPopulationSize',200,@(x)x>=1);
parser.addParameter('gaIterations',500,@(x)x>=1);
parser.addParameter('fminconIterations',500,@(x)x>=1);
parser.addParameter('verbose',true,@islogical);
parser.addParameter('weighting',[],@isstruct);
parser.parse(varargin{:});
arg = parser.Results;
initial = parser.Results.initial;
eps = parser.Results.eps;
lb = parser.Results.lb;
ub = parser.Results.ub;
fix = parser.Results.fix;
Usep = parser.Results.Usep;
w = parser.Results.w;
popSizeGA = parser.Results.gaPopulationSize;
iterGA = parser.Results.gaIterations;
iterFmincon = parser.Results.fminconIterations;
verbose = parser.Results.verbose;
weighting = parser.Results.weighting;
if ~isstruct(eps)
    e = cell(1,length(msmrFieldnames));
    for k = 1:length(e); e{k} = eps; end
    eps = cell2struct(e,msmrFieldnames,2);
end

% Calculate the number of degrees of freedom in the optimization
% (total number of MSMR variables to optimize).
ndim = sum([pack(msmrFields,fix).len]);

% Define optimization bounds. Ensure all parameters >0.
for k = 1:length(msmrFieldnames)
    field = msmrFieldnames{k};
    if isfield(fix,field)
        % No need to define upper and lower bounds for
        % fixed (non-optimized) parameters.
        continue;
    end
    if isfield(lb,field)
        lb.(field) = ones(msmrFields.(field).len,1).*lb.(field);
    elseif ~isempty(initial)
        lb.(field) = initial.(field).*(1-eps.(field));
    else
        lb.(field) = zeros(msmrFields.(field).len,1);
    end
    if isfield(ub,field)
        ub.(field) = ones(msmrFields.(field).len,1).*ub.(field);
    elseif ~isempty(initial)
        ub.(field) = initial.(field).*(1+eps.(field));
    else
        ub.(field) = Inf*ones(msmrFields.(field).len,1);
    end
    lb.(field) = max(0,lb.(field));
    ub.(field) = max(0,ub.(field));
end

% Ensure compositions remain between 0 and 1.
for k = 1:length(msmrFieldnamesComposition)
    field = msmrFieldnamesComposition{k};
    if isfield(fix,field)
        % No need to define upper and lower bounds for
        % fixed (non-optimized) parameters.
        continue;
    end
    lb.(field) = max(min(lb.(field),1),0);
    ub.(field) = max(min(ub.(field),1),0);
end

% Specify initial population for the genetic algorithm.
population = zeros(popSizeGA,ndim);
for kmbr = 1:popSizeGA
    % Gallery potentials.
    if ~isfield(fix,'U0')
        mbr.U0 = lb.U0 + rand(J,1).*(ub.U0-lb.U0);
        mbr.U0 = sort(mbr.U0,'descend'); % ensure U1 > U2 > U3 > ...
    end
    % Gallery partial compositions.
    if ~isfield(fix,'X')
        mbr.X = lb.X + rand(J,1).*(ub.X-lb.X);
        mbr.X = mbr.X/sum(mbr.X);  % ensure sum is 1
    end
    % Gallery disorder factors.
    if ~isfield(fix,'omega')
        mbr.omega = lb.omega + rand(J,1).*(ub.omega-lb.omega);
    end
    % Minumum/maximum composition.
    if ~isfield(fix,'thetamin')
        mbr.thetamin = lb.thetamin + rand(1,1)*(ub.thetamin-lb.thetamin);
    end
    if ~isfield(fix,'thetamax')
        mbr.thetamax = lb.thetamax + rand(1,1)*(ub.thetamax-lb.thetamax);
    end
    population(kmbr,:) = pack(mbr,fix).';
end % for
if ~isempty(initial)
    % Add initial value to population.
    [~,ind] = sort(initial.U0,'descend');  % ensure U1 > U2 > U3 > ...
    initial.U0 = initial.U0(ind);
    initial.X = initial.X(ind);
    initial.omega = initial.omega(ind);
    population(1,:) = pack(initial,fix).';
end

% Convert bounds to vectors.
lb = pack(lb,fix); 
ub = pack(ub,fix);

% Define the U1 > U2 > U3 > ... constraint. Ensures the
% galleries do not "swap" positions so that ub / lb constraints may be
% applied to each individual gallery.
orderU.U0 = toeplitz([-1 zeros(1,J-2)],[-1 1 zeros(1,J-2)])';
fnames = setdiff(msmrFieldnames,{'U0'});
for k = 1:length(fnames)
    field = fnames{k};
    orderU.(field) = zeros(msmrFields.(field).len,J-1);
end
A = pack(orderU,fix)';
B = -ones(J-1,1).*Usep(:);  % Enforce minimum gallery potential separation.

% Define the sum(Xj)=1 constraint.
sumX.X = ones(J,1);
fnames = setdiff(msmrFieldnames,{'X'});
for k = 1:length(fnames)
    field = fnames{k};
    sumX.(field) = zeros(msmrFields.(field).len,1);
end
Aeq = pack(sumX,fix)';
Beq = 1;

% Configure and perform hybrid Genetic Algorithm optimization.
if verbose
    display = 'iter';
else
    display = 'off';
end
optionsHybrid = optimoptions(@fmincon, ...
    'Display', display,...
    'MaxFunEvals', 1e6, ...
    'MaxIter', iterFmincon,...
    'TolFun', 1e-10, ...
    'TolX', 1e-10, ...
    'TolCon', 1e-10);
options = optimoptions(@ga,...
    'Display', display, ...
    'PopulationSize', popSizeGA, ...
    'InitialPopulationMatrix', population, ...
    'MaxGenerations', iterGA, ...
    'FunctionTolerance', 1e-10,...
    'MaxStallGenerations', 100,...
    'HybridFcn', {@fmincon,optionsHybrid});
[vect, Jcost] = ga(@cost,ndim,A,B,Aeq,Beq,lb,ub,[],options);
est = unpack(vect,fix);

% Sort galleries *ascending* by U0.
[~,ind] = sort(est.U0);
est.U0    = est.U0(ind);
est.X     = est.X(ind);
est.omega = est.omega(ind);

% Collect output.
data.est = est;
data.Jcost = Jcost;
data.origin__ = 'fitMSMR';
data.arg__ = arg;

function Jcost = cost(vect)
    % Unpack parameter vector into MSMR model.
    params = unpack(vect,fix);

    % Accumulate cost as differences between model prediction
    % and laboratory measurements (more than one
    % laboratory-derived OCP estimate may be supplied).
    Jcost = 0;
    for idxEstimate = 1:length(ocp)
        TdegC = ocp(idxEstimate).TdegC;
        Ulab = ocp(idxEstimate).U;
        Zlab = ocp(idxEstimate).Z;
        dZlab = abs(ocp(idxEstimate).dZ); % !!! use absolute diff cap

        % Evalulate MSMR model over same Uocp vector as lab data.
        % (We have to re-evalulate the model prediction for each
        % OCP estimate, as the temperature T may differ.)
        U0 = params.U0;
        X = params.X;
        omega = params.omega;
        Umod = Ulab(:)';
        f = F/R/(TdegC+273.15);
        gj = exp(f*(Umod-U0)./omega);
        xj = X./(1+gj);
        Zmod = sum(xj,1);
        dZmod = f*sum((X./omega).*gj./(1+gj).^2,1); % !!! use absolute diff cap

        % Convert model vectors over absolute composition to relative
        % composition for comparison with data.
        diffzlab = max(Zlab) - min(Zlab);
        Zmod = (Zmod - params.thetamin)/(params.thetamax - params.thetamin);
        Zmod = min(Zlab) + Zmod*diffzlab;
        dZmod = dZmod*diffzlab/(params.thetamax - params.thetamin);

        % Compute residuals.
        zResid = Zmod(:) - Zlab(:);
        dzResid = dZmod(:) - dZlab(:);

        % Apply weighting.
        for weight = weighting(:)'
            idx = weight.getInterval(Ulab);
            zResid(idx) = zResid(idx)*weight.multiplier;
            dzResid(idx) = dzResid(idx)*weight.multiplier;
        end

        % Compute cost.
        Jcost = Jcost + sum(zResid.^2) + w*sum(dzResid.^2);
    end
end

end


% Utility functions -------------------------------------------------------

function vect = pack(params, fixed)
%PACK Reduce an MSMR parameter structure to a vector of parameters. 
%
% VECT = PACK(MODEL) compacts parameters in MODEL into a 
%   column vector. MODEL may be an instance of MSMR or a
%   structure containing fields U0,X,omega,thetamin,thetamax.
%
% VECT = PACK(...,FIXEDPARAMS) excludes the fields of the
%   structure FIXEDPARAMS from the parameter vector.

if ~exist('fixed','var')
    fixed = struct;
end

if isfield(fixed,'X'), X = []; else, X = params.X; end
if isfield(fixed,'U0'), U0 = []; else, U0 = params.U0; end
if isfield(fixed,'omega'), omega = []; else, omega = params.omega; end
if isfield(fixed,'thetamin'), thetamin = []; else, thetamin = params.thetamin; end
if isfield(fixed,'thetamax'), thetamax = []; else, thetamax = params.thetamax; end

vect = [X; U0; omega; thetamin; thetamax];
end

function params = unpack(vect, fixed)
%UNPACK Restore a MSMR parameter structure from the given vector.

if ~exist('fixed','var')
    fixed = struct;
end

numVects = 3; numScalars = 2;
if isfield(fixed,'X'), numVects = numVects - 1; X = fixed.X; J = length(fixed.X); end
if isfield(fixed,'U0'), numVects = numVects - 1; U0 = fixed.U0; J = length(fixed.U0); end
if isfield(fixed,'omega'), numVects = numVects - 1; omega = fixed.omega; J = length(fixed.omega); end
if isfield(fixed,'thetamin'), numScalars = numScalars - 1; thetamin = fixed.thetamin; end
if isfield(fixed,'thetamax'), numScalars = numScalars - 1; thetamax = fixed.thetamax; end

cursor = 1;
if ~exist('J','var'), J = (length(vect)-numScalars)/numVects; end
if ~exist('X','var'), X = vect(cursor:J+cursor-1); cursor = cursor + J; end
if ~exist('U0','var'), U0 = vect(cursor:J+cursor-1); cursor = cursor + J; end
if ~exist('omega','var'), omega = vect(cursor:J+cursor-1); cursor = cursor + J; end
if ~exist('thetamin','var'), thetamin = vect(cursor); cursor = cursor + 1; end
if ~exist('thetamax','var'), thetamax = vect(cursor); end

% Reconstruct MSMR model.
params.U0 = U0(:);
params.X = X(:);
params.omega = omega(:);
params.thetamin = thetamin;
params.thetamax = thetamax;
end