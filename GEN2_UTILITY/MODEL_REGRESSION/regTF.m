function regspec = regTF(regtype,eis,ocp,init,varargin)
%REGTF Generate specification (bounds, initial value, etc.) for TF model 
% regression.
%
% -- Usage --
% regspec = REGTF(regtype,eis,ocp,init)
% regspec = REGTF(...,'fix',fix)
% regspec = REGTF(...,'TrefdegC',Tref)
%
% -- Input --
% regtype = regression type, choose from the following:
%   'msmr' use MSMR models for Ds and i0
%   'lut'  use lookup table models for Ds and i0
%
% eis = struct or struct array (one entry for each temperature) of the
% following fields:
%   .socPct    SOC vector (true SOC computed by coulomb-counting) [%]
%   .freq      Frequency vector [Hz]
%   .Z         Linear impedance matrix (dim1<=>freq, dim2<=>socPct) [Ohm]
%   .TdegC     Temperature [degC]
%
% ocp = struct with the following fields:
%   .U0        Column vector of MSMR { U0_j  } parameter values
%   .X         Column vector of MSMR { X_j   } parameter values
%   .omega     Column vector of MSMR { omega } parameter values
%   .theta0    Porous electrode composition at 0% SOC
%   .theta100  Porous electrode composition at 100% SOC
%   .Q         Total capacity [Ah]
%
% init = cell model of initial parameter values (initial guess); obtain
%   from loadCellModel.m
%
% fix = struct of model parameters to fix to certain values (overrides
%   defaults)
%
% Tref = reference temperature to use in Arrhenius relations for
%   temperature-dependent parameters [degC]

parser = inputParser;
parser.addRequired('regtype',@(x)any(strcmpi(x,{'msmr','lut'})));
parser.addRequired('eis',@(x)isvector(x)&&isstruct(x));
parser.addRequired('ocp',@(x)isscalar(x)&&isstruct(x));
parser.addRequired('init',@(x)isstruct(x)&&isscalar(x));
parser.addParameter('TrefdegC',25,@(x)isscalar(x)&&isnumeric(x));
parser.addParameter('fix',[],@(x)isstruct(x)&&isscalar(x)||isempty(x));
parser.parse(regtype,eis,ocp,init,varargin{:});
arg = parser.Results; % struct of validated arguments

% Constants.
tplus0 = 0.4;    % Li+ transference number, guess for estimating psi
R = TB.const.R;  % molar gas constant [J/mol/K]
F = TB.const.F;  % Faraday's constant [C/mol]

% Determine the multiplicity of the dataset, i.e. the number of
% temperatures included.
ntemp = length(arg.eis);

% Determine number of SOC setpoints at each temperature, minimum SOC, and
% maximum SOC.
nsoc = arrayfun(@(e)length(e.socPct),eis);
nsocavg = round(mean(nsoc));
minsocPct = min(arrayfun(@(e)min(e.socPct),eis));
maxsocPct = max(arrayfun(@(e)min(e.socPct),eis));

% Determine min/max composition of porous electrode.
theta1 = arg.ocp.theta0 + (minsocPct/100)*(arg.ocp.theta100-arg.ocp.theta0);
theta2 = arg.ocp.theta0 + (maxsocPct/100)*(arg.ocp.theta100-arg.ocp.theta0);
mintheta = min(theta1,theta2);
maxtheta = max(theta1,theta2);

% Generate composition vector for LUT for Ds and i0.
thetaLUT = linspace(mintheta,maxtheta,nsocavg).';

% Determine number of MSMR galleries.
J = length(arg.ocp.U0);

% Form temperature vector.
TdegC = [arg.eis.TdegC].';


% Build model -------------------------------------------------------------
% Known parameters.
params.const.Q = fastopt.param('fix',arg.ocp.Q);
params.pos.theta0 = fastopt.param('fix',arg.ocp.theta0);
params.pos.theta100 = fastopt.param('fix',arg.ocp.theta100);
params.pos.X = fastopt.param('fix',arg.ocp.X);
params.pos.U0 = fastopt.param('fix',arg.ocp.U0);
params.pos.omega = fastopt.param('fix',arg.ocp.omega);

% Fixed parameters.
params.neg.Rf = fastopt.param('fix',0);  % Lump into tab resistance.
% Symmetry factors for the positive electrode are not very identifiable,
% we'll asume alpha=0.5 for all galleries.
if strcmpi(arg.regtype,'lut')
    params.pos.alphaLinear = fastopt.param('fix',0.5*ones(length(thetaLUT),1));
else
    params.pos.alpha = fastopt.param('fix',0.5*ones(J,1));
end
% Symmetry factor for the negative electrode is not identifiable due to
% linear nature of the small-signal impedance model.
params.neg.alpha = fastopt.param('fix',0.5);
% Fix composition of interpolation points for k0 and Ds.
if strcmpi(arg.regtype,'lut')
    params.pos.k0Theta = fastopt.param('fix',thetaLUT);
    params.pos.DsTheta = fastopt.param('fix',thetaLUT);
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
if strcmpi(arg.regtype,'lut')
    params.pos.DsLinear = fastopt.param( ...
        'len',length(thetaLUT),'logscale',true,'tempfcn','Eact');
    params.pos.k0Linear = fastopt.param( ...
        'len',length(thetaLUT),'logscale',true,'tempfcn','Eact');
else
    params.pos.Dsref = fastopt.param( ...
        'logscale',true,'tempfcn','Eact');
    params.pos.k0 = fastopt.param( ...
        'len',J,'logscale',true,'tempfcn','Eact');
end
params.pos.nF = fastopt.param;
params.pos.tauF = fastopt.param('logscale',true);
params.pos.Cdl = fastopt.param;
params.pos.nDL = fastopt.param;
params.pos.tauDL = fastopt.param('logscale',true);
params.pos.Rf = fastopt.param;
params.pos.Rdl = fastopt.param;
params.pos.sigma = fastopt.param;

% Lithium-metal electrode parameters.
params.neg.k0 = fastopt.param('logscale',true,'tempfcn','Eact');
params.neg.Cdl = fastopt.param;
params.neg.nDL = fastopt.param;
params.neg.tauDL = fastopt.param('logscale',true);
params.neg.Rdl = fastopt.param;

% Cell package parameters.
params.pkg.R0 = fastopt.param('tempfcn','lut');  % Tab resistance.
params.pkg.L0 = fastopt.param('tempfcn','lut');  % Package/cable inductace.

% Load additional fixed parameters.
% (Overwrite defaults above).
if ~isempty(arg.fix)
    regs = fieldnames(arg.fix);
    for kr = 1:length(regs)
        reg = regs{kr};
        paramnames = fieldnames(arg.fix.(reg));
        for kp = 1:length(paramnames)
            pname = paramnames{kp};
            params.(reg).(pname) = fastopt.param('fix',arg.fix.(reg).(pname));
        end
    end
end

modelspec = fastopt.modelspec(params, ...
    'tempsdegC',TdegC,'TrefdegC',arg.TrefdegC);


% Build initial guess -----------------------------------------------------

% Convert cell model to reduced-layer Warburg-resistance model.
% Also convert the kinetics and diffusion sub-models.
if strcmpi(arg.regtype,'lut')
    kStruct.type = 'linear';
    kStruct.theta = thetaLUT;  % data point near each lab SOC setpoint
else
    kStruct.type = 'msmr';
end
if strcmpi(arg.regtype,'lut')
    dStruct.type = 'linear';
    dStruct.theta = thetaLUT;  % data point near each lab SOC setpoint
else
    dStruct.type = 'msmr';
end
initialModel = convertCellModel(arg.init,'RLWRM', ...
    'KineticsModel',kStruct,'SolidDiffusionModel',dStruct);

% Load OCP parameters into initial model.
p = struct;
p.const.Q = arg.ocp.Q;
p.pos = rmfield(arg.ocp,'Q');
initialModel = setCellParam(initialModel,p);

% Fetch initial values of model parameters.
initial = getCellParams(initialModel,'TdegC',arg.TrefdegC);

% !!! Important: Lump Rf(neg) with R0(pkg)!
initial.pkg.R0 = initial.pkg.R0 + initial.neg.Rf;

% Generate initial guess; pack/unpack to set values of fixed parameters.
% Add activation energy parameters.
initial.pos.tauW_Eact = 1000;
initial.pos.kappa_Eact = 1000;
initial.eff.tauW_Eact = 1000;
initial.eff.kappa_Eact = 1000;
if strcmpi(arg.regtype,'lut')
    initial.pos.DsLinear_Eact = 1000;
    initial.pos.k0Linear_Eact = 1000;
else
    initial.pos.Dsref_Eact = 1000;
    initial.pos.k0_Eact = 1000;
end
initial.neg.k0_Eact = 1000;
init = fastopt.unpack(fastopt.pack(initial,modelspec),modelspec);


% Build lower/upper bounds ------------------------------------------------

% Electrolyte parameters.
lb.const.W = 0.1;                   ub.const.W = 20;
lb.eff.tauW = init.eff.tauW/100;    ub.eff.tauW = init.eff.tauW*100;
lb.eff.tauW_Eact = 0;               ub.eff.tauW_Eact = 100000;
lb.pos.tauW = init.pos.tauW/100;    ub.pos.tauW = init.pos.tauW*100;
lb.pos.tauW_Eact = 0;               ub.pos.tauW_Eact = 100000;
lb.eff.kappa = init.eff.kappa/100;  ub.eff.kappa = init.eff.kappa*100;
lb.eff.kappa_Eact = 0;              ub.eff.kappa_Eact = 100000;
lb.pos.kappa = init.pos.kappa/100;  ub.pos.kappa = init.pos.kappa*100;
lb.pos.kappa_Eact = 0;              ub.pos.kappa_Eact = 100000;

% Porous-electrode parameters.
if strcmpi(arg.regtype,'msmr')
    lb.pos.Dsref = 1e-9;            ub.pos.Dsref = 1e3;
    lb.pos.Dsref_Eact = 0;          ub.pos.Dsref_Eact = 100000;
else
    lb.pos.DsLinear = 1e-5*ones(length(thetaLUT),1);  % lower bound
    ub.pos.DsLinear = 1*ones(length(thetaLUT),1);     % upper bound
    lb.pos.DsLinear_Eact = 0;       ub.pos.DsLinear_Eact = 100000;
end
lb.pos.nF = 0.5;                    ub.pos.nF = 1;
lb.pos.tauF = 0.1;                  ub.pos.tauF = 10000;
if strcmpi(arg.regtype,'msmr')
    lb.pos.k0 = 1e-4*ones(J,1);     ub.pos.k0 = 100*ones(J,1);
    lb.pos.k0_Eact = 0;             ub.pos.k0_Eact = 100000;
else
    lb.pos.k0Linear = 1e-4*ones(length(thetaLUT),1);  % lower bound
    ub.pos.k0Linear = 100*ones(length(thetaLUT),1);   % upper bound
    lb.pos.k0Linear_Eact = 0;       ub.pos.k0Linear_Eact = 100000;
end
lb.pos.sigma = init.pos.sigma/1000; ub.pos.sigma = init.pos.sigma*1000;
lb.pos.Cdl = init.pos.Cdl/10;       ub.pos.Cdl = init.pos.Cdl*100;
lb.pos.nDL = 0.7;                   ub.pos.nDL = 1;
lb.pos.tauDL = 1e-6;                ub.pos.tauDL = 1e6;
lb.pos.Rf = init.pos.Rf/100;        ub.pos.Rf = init.pos.Rf*100;
lb.pos.Rdl = init.pos.Rdl/100;      ub.pos.Rdl = init.pos.Rdl*100;

% Lithium-metal electrode parameters.
lb.neg.k0 = init.neg.k0/1000;       ub.neg.k0 = init.neg.k0*1000;
lb.neg.k0_Eact = 0;                 ub.neg.k0_Eact = 100000;
lb.neg.Cdl = init.neg.Cdl/10;       ub.neg.Cdl = init.neg.Cdl*10;
lb.neg.nDL = 0.7;                   ub.neg.nDL = 1;
lb.neg.tauDL = 1e-6;                ub.neg.tauDL = 1e6;
lb.neg.Rdl = init.neg.Rdl/100;      ub.neg.Rdl = init.neg.Rdl*100;

% Cell package parameters.
lb.pkg.R0 = 0;                      ub.pkg.R0 = init.pkg.R0*10;
lb.pkg.L0 = 0;                      ub.pkg.L0 = init.pkg.L0*100;


% Collect output ----------------------------------------------------------
regspec.type = arg.regtype;
regspec.ntemp = ntemp;
regspec.TdegC = TdegC;
regspec.eis = arg.eis;
regspec.ocp = arg.ocp;
regspec.modelspec = modelspec;
regspec.init = init;
regspec.lb = lb;
regspec.ub = ub;
regspec.origin__ = 'regTF';
regspec.arg__ = arg;

end