function data = tfLMB(S, model, varargin)
%TFLMB Compute the transfer-functions for a full LMB cell. Linear and
%  second-harmonic response. Faster than using tfXX functions for computing
%  cell impedance.
%
% data = TFLMB(S,model) constructs functions for evalulating the LMB
%   transfer functions at the frequency points specified by the vector S.
%   The output DATA is a structure with the following fields:
%      tfThetae(X): a function that returns the values of the
%                   Thetae(s)/Iapp(s) transfer function at the positions
%                   specified in the vector X. In particular, the output is
%                   an ns-by-nx matrix where nx is the number of positions
%                   in the vector X.
%
% data = TFLMB(...,'TdegC',TdegC) performs the computation at temperature 
%    TdegC instead of the default 25degC.
%
% data = TFLMB(...,'ConsolidateLayers','seponly') ignores the dead-lithium
%   layer, taking into account the separator parameters only.
%
% data = TFLMB(...,'ConsolidateLayers','combine') combines the dead-lithium
%   and separator layers into a single chemically-intert porous layer.
%
% data = TFLMB(...,'Calc11',false) disables calculation of the linear
%   response; only the parameters evalulated at the setpoint are returned.
%
% data = TFLMB(...,'Calc22',true) also calculates the second-harmonic
%   response. The linear response will also be computed, as it is needed to
%   find the second-harmonic response 

parser = inputParser;
parser.StructExpand = false;
parser.addRequired('S',@(x)isvector(x));
parser.addRequired('model',@(x)isstruct(x));
parser.addParameter('TdegC',25,@isscalar);
parser.addParameter('socPct',50,@isscalar);
parser.addParameter('Calc11',true,@islogical);
parser.addParameter('Calc22',false,@islogical);
parser.addParameter('ParameterValues',[],@isstruct);
parser.parse(S,model,varargin{:});
TdegC = parser.Results.TdegC;
socPct = parser.Results.socPct;
calc11 = parser.Results.Calc11;
calc22 = parser.Results.Calc22;
param = parser.Results.ParameterValues;
S = S(:);  % force S to be a column vector

if isempty(param)
    param = getParameterValues(model,TdegC,socPct,S);
end
if calc11 || calc22
    data.h11 = tfLMB11(param,S,param.layerReduction);
end
if calc22
    data.h22 = tfLMB22(data.h11,param.layerReduction);
end
data.param = param;

end % tfLMB()

function p = getParameterValues(model,TdegC,socPct,S)
%GETPARAMETERVALUES

T = TdegC+273.15;
R = 8.3144598;      % Molar gas constant [J/mol K]
F = 96485.3329;     % Faraday constant [C/mol]

% Get lithiation of positive electrode.
if isfield(model,'function')
    % Toolbox cell model.
    theta0p = model.function.pos.theta0();
    theta100p = model.function.pos.theta100();
else
    % Set of model parameters already evalulated at setpoint.
    theta0p = model.pos.theta0;
    theta100p = model.pos.theta100;
end
p.socp = socPct/100;
p.thetap = theta0p + p.socp*(theta100p - theta0p);
thetap = p.thetap;

if isfield(model,'function')
    % Toolbox cell model.
    getReg = @(reg)model.function.(reg);
    isReg = @(reg)isfield(model.function,reg);
    getParam = @(reg,p)model.function.(reg).(p)(thetap,T);
    isParam = @(reg,p)isfield(model.function.(reg),p);
else
    % Set of model parameters already evalulated at setpoint.
    getReg = @(reg)model.(reg);
    isReg = @(reg)isfield(model,reg);
    getParam = @(reg,p)model.(reg).(p);
    isParam = @(reg,p)isfield(model.(reg),p);
end

% Collect parameters.
p.T           = T;
p.TdegC       = TdegC;
p.R           = R;
p.F           = F;
p.psi         = getParam('const','psi');
p.W           = getParam('const','W');
p.Q           = getParam('const','Q');  
p.sigmap      = getParam('pos','sigma');
p.kappap      = getParam('pos','kappa');
p.taup        = getParam('pos','tauW');
p.nF          = getParam('pos','nF');
p.Rfp         = getParam('pos','Rf');
p.k0p         = getParam('pos','k0');
p.alphap      = getParam('pos','alpha');
p.Rdlp        = getParam('pos','Rdl');
p.Cdlp        = getParam('pos','Cdl');
p.nDLp        = getParam('pos','nDL');
p.theta0p     = getParam('pos','theta0');
p.theta100p   = getParam('pos','theta100');
p.Rfn         = getParam('neg','Rf');
p.Rdln        = getParam('neg','Rdl');
p.Cdln        = getParam('neg','Cdl');
p.nDLn        = getParam('neg','nDL');
p.k0n         = getParam('neg','k0');
p.alphan      = getParam('neg','alpha');
% Code below implements Warburg CPE factor for inert layers.
p.nEs         = 1;
p.nEdl        = 1;

% Determine which electrolyte layers are present; extract parameters.
layerReduction = 'off';
if isReg('eff')
    p.kappas    = getParam('eff','kappa');
    p.taus      = getParam('eff','tauW');
    p.kappad    = Inf;
    p.taud      = 0;
    layerReduction = 'combine';
elseif isReg('dll') && isReg('sep')
    p.kappas    = getParam('sep','kappa');
    p.taus      = getParam('sep','tauW');
    p.kappad    = getParam('dll','kappa');
    p.taud      = getParam('dll','tauW');
else
    error(['Cell model must inlcude dll and sep regions ' ...
        'or single eff region.']);
end
% Legcacy code expects qe, kD electrolyte parameters; calculate them
% from tau, W.
p.qep = R*p.kappap*p.taup/3600/F;
p.qes = R*p.kappas*p.taus/3600/F;
p.qed = R*p.kappad*p.taud/3600/F;
p.kD = -p.psi*p.W;
if isnan(p.qed)
    p.qed = 0; % dll combined with sep
end

% Get Uocp, dUocp, and Rct for positive electrode.
if isParam('pos','Uocp') && isParam('pos','dUocp')
    % Non-MSMR model.
    p.Uocpp = getParam('pos','Uocp');
    p.dUocpp = getParam('pos','dUocp');
    if isParam('pos','d2Uocp')
        p.d2Uocpp = getParam('pos','d2Uocp');
    elseif calc22
        if isfield(model,'function')
            % We need d2Uocp to compute the second-harmonic impedance.
            % Use first-difference approximation.
            warning(['d2Uocp(pos) is missing from cell model, but it is required ' ...
                'to evalulate the second-harmonic impedance. Using ' ...
                'first-difference approximation of dUocp(pos).']);
            thetaVect = linspace(p.theta100p,p.theta0p,10000);
            dUocpVect = model.function.pos.dUocp(thetaVect,T);
            d2thetaVect =  thetaVect(1:end-1) + diff(thetaVect)/2;
            d2UocpVect = diff(dUocpVect)./diff(thetaVect);
            p.d2Uocpp = interp1(d2thetaVect,d2UocpVect,p.thetap,'linear','extrap');
        else
            % dUocp function not available to perform first-difference
            % aproximation of d2Uocp.
            error(['d2Uocp(pos) is missing from cell model, but is required ' ...
                'to compute the second-harmonic impedance!']);
        end
    end
    p.i0p = p.k0p.*(1-p.thetap).^(1-p.alphap).*p.thetap.^(p.alphap);
    p.Rctp = R*T/p.i0p/F;
    p.Rct2invp = F*p.i0p*((1-p.alphap).^2-p.alphap.^2)/R/T;
elseif isParam('pos','X') && isParam('pos','U0') && isParam('pos','omega')
    % MSMR model.
    electrode = MSMR(getReg('pos'),'sortParams',false);
    [RctVect, UocpVect, dUocpVect, thetaVect, dataRct] = electrode.Rct( ...
        getReg('pos'),'npoints',10000,'TdegC',TdegC);
    alphap = p.alphap(:);
    Rct2invVect = sum(dataRct.i0j.*((1-alphap).^2-alphap.^2))*dataRct.f;
    p.Uocpp = interp1(thetaVect,UocpVect,p.thetap,'linear','extrap');
    p.dUocpp = interp1(thetaVect,dUocpVect,p.thetap,'linear','extrap');
    p.d2Uocpp = interp1(thetaVect,dataRct.d2Uocp,p.thetap,'linear','extrap');
    p.Rctp = interp1(thetaVect,RctVect,p.thetap,'linear','extrap');
    p.Rct2invp = interp1(thetaVect,Rct2invVect,p.thetap,'linear','extrap');
else
    error(['Neither the Uocp(pos) and dUocp(pos) functions ' ...
        'or MSMR parameters are provided in the cell model. ' ...
        'Please supply one or the other.']);
end

% Calculate Rct for negative electrode.
% (Assumes thetas=1 and thetae=1.)
p.Rctn = R*T/p.k0n/F;
p.Rct2invn = F*p.k0n*((1-p.alphan).^2-p.alphan.^2)/R/T;

% Get solid diffusivity.
if isParam('pos','Ds')
    p.Dsp = getParam('pos','Ds');
elseif isParam('pos','Dsref')
    Dsrefp = getParam('pos','Dsref');
    p.Dsp = -Dsrefp*F/R/T*p.thetap*(1-p.thetap)*p.dUocpp;
else
    error(['Neither Ds nor Dsref provided by the cell model. ' ...
        'Please supply one or the other.']);
end

% Compute Zse.
DeltaQp = abs(p.theta100p - p.theta0p);
p.csmaxp = 10800*p.Q*p.Dsp/DeltaQp;
p.gammap = sqrt((S/p.Dsp).^p.nF);
p.Zsp = (p.dUocpp/p.csmaxp)...
   *((p.gammap.^2 + 3*(1-p.gammap.*coth(p.gammap)))...
   ./(p.gammap.^2.*(1-p.gammap.*coth(p.gammap))) - 3*p.Dsp./S);
p.Zdlp = p.Rdlp + 1./p.Cdlp./(S.^p.nDLp);
p.Zsep = p.Rfp + 1./(1./(p.Rctp + p.Zsp) + 1./p.Zdlp);
p.Zdln = p.Rdln + 1./p.Cdln./(S.^p.nDLn);
p.Zsen = p.Rfn + p.Zdln.*p.Rctn./(p.Zdln+p.Rctn);

% Compute Zse2 (double the frequency!)
p.gamma2p = sqrt((2*S/p.Dsp).^p.nF);
p.Zs2p = (p.dUocpp/p.csmaxp)...
   *((p.gamma2p.^2 + 3*(1-p.gamma2p.*coth(p.gamma2p)))...
   ./(p.gamma2p.^2.*(1-p.gamma2p.*coth(p.gamma2p))) - 3*p.Dsp./S/2);
p.Zdl2p = p.Rdlp + 1./p.Cdlp./(2*S).^p.nDLp;
p.Zse2p = p.Rfp + p.Zdl2p.*(p.Rctp+p.Zs2p)./(p.Zdl2p+p.Rctp+p.Zs2p);
p.Zdl2n = p.Rdln + 1./p.Cdln./(2*S).^p.nDLn;
p.Zse2n = p.Rfn + p.Zdl2n.*p.Rctn./(p.Zdl2n+p.Rctn);

% Additional.
p.split = splitFactory();
p.xlim = xlimFactory();
p.layerReduction = layerReduction;

    function fcn = splitFactory()
        function [Z, bins] = split1(X)
            X = X(:).';
            bins.n = X==0;
            bins.d = 0<X&X<=1;
            bins.s = 1<X&X<2;
            bins.p = 2<=X&X<=3;  % Some quantities not defined in s, so ensure x=2 is in p!
            Z.n = X(bins.n);
            Z.d = X(bins.d);
            Z.s = X(bins.s)-1;
            Z.p = 3-X(bins.p);
        end
        function [Z, bins] = split2(X)
            X = X(:).';
            bins.n = X==0;
            bins.d = false(size(X));
            bins.s = 0<X&X<1;
            bins.p = 1<=X&X<=2;  % Some quantities not defined in s, so ensure x=1 is in p!
            Z.n = X(bins.n);
            Z.d = [];
            Z.s = X(bins.s);
            Z.p = 2-X(bins.p);
        end
        if strcmpi(layerReduction,'off')
            % Two inert layers (dll and sep).
            fcn = @split1;
        else
            % One inert layer (sep).
            fcn = @split2;
        end
    end
    function x = xlimFactory()
        if strcmpi(layerReduction,'off')
            % Two inert layers (dll and sep).
            x.min = 0;
            x.sp = 2;
            x.max = 3;
        else
            % One inert layer (sep).
            x.min = 0;
            x.sp = 1;
            x.max = 2;
        end
    end
end % getParameterValues()

% First-Harmonic ----------------------------------------------------------

function data = tfLMB11(p,S,layerReduction)
    ns = length(S);

    % Compute quantities needed to solve for electrolyte concentration
    % coefficients.
    mu1p = (1/p.kappap/(p.psi*p.T))./p.Zsep;
    mu2p = (p.kD*p.T/(p.psi*p.T)/p.kappap)./p.Zsep - ((3600*p.qep/(p.psi*p.T)/p.kappap).*S);
    tau1p = (1/p.sigmap + 1/p.kappap)./p.Zsep - mu2p; 
    tau2p = ((3600*p.qep/(p.psi*p.T)/p.kappap).*S).*(1/p.sigmap + 1/p.kappap)./p.Zsep; 
    Lambda1P = sqrt(0.5*(tau1p-sqrt(tau1p.^2-4*tau2p))); 
    Lambda2P = sqrt(0.5*(tau1p+sqrt(tau1p.^2-4*tau2p))); 
    lambda1p = Lambda1P.^3+mu2p.*Lambda1P;
    lambda2p = Lambda2P.^3+mu2p.*Lambda2P;
    expnLambda1P = exp(-Lambda1P);
    expnLambda2P = exp(-Lambda2P);
    LambdaS = (3600*p.qes*S/p.psi/p.T/p.kappas).^(p.nEs/2);
    expnLambdaS = exp(-LambdaS);
    if strcmpi(layerReduction,'off')
        LambdaD = (3600*p.qed*S/p.psi/p.T/p.kappad).^(p.nEdl/2);
        expnLambdaD = exp(-LambdaD);
    else
        LambdaD = [];
        expnLambdaD = [];
    end
    
    meta.aux.mu1p = mu1p;
    meta.aux.mu2p = mu2p;
    meta.aux.tau1p = tau1p;
    meta.aux.tau2p = tau2p;
    
    meta.S = S;
    meta.ns = ns;
    meta.TdegC = p.TdegC;
    meta.T = p.T;
    meta.R = p.R;
    meta.F = p.F;
    meta.param = p;
    meta.LambdaS = LambdaS;
    meta.LambdaD = LambdaD;
    meta.Lambda1P = Lambda1P;
    meta.Lambda2P = Lambda2P;
    
    if strcmpi(layerReduction,'off')
        C = solveTwoInertLayers;
        meta.c1d = C(:,1);
        meta.c2d = C(:,2);
        meta.c1s = C(:,3);
        meta.c2s = C(:,4);
        meta.c1p = C(:,5);
        meta.c2p = C(:,6);
        meta.c3p = C(:,7);
        meta.c4p = C(:,8);
    else
        C = solveSingleInertLayer;
        meta.c1s = C(:,1);
        meta.c2s = C(:,2);
        meta.c1p = C(:,3);
        meta.c2p = C(:,4);
        meta.c3p = C(:,5);
        meta.c4p = C(:,6);
    end
    
    meta.j1p = meta.c1p.*(-(p.psi*p.T)*p.kappap*meta.Lambda1P.^2+3600*p.qep.*S);
    meta.j2p = meta.c2p.*(-(p.psi*p.T)*p.kappap*meta.Lambda1P.^2+3600*p.qep.*S);
    meta.j3p = meta.c3p.*(-(p.psi*p.T)*p.kappap*meta.Lambda2P.^2+3600*p.qep.*S);
    meta.j4p = meta.c4p.*(-(p.psi*p.T)*p.kappap*meta.Lambda2P.^2+3600*p.qep.*S);
    meta.divFactorIfp = p.Zdlp./(p.Rctp+p.Zsp+p.Zdlp);
    meta.divFactorIdlp = (p.Rctp+p.Zsp)./(p.Rctp+p.Zsp+p.Zdlp);
    meta.divFactorIfn = p.Zdln./(p.Rctn+p.Zdln);
    meta.divFactorIdln = p.Rctn./(p.Rctn+p.Zdln);
    
    data.tfThetae = @(X)getThetae(meta,X,0);
    data.tfDThetae = @(r,X)getThetae(meta,X,r);
    data.tfIfdl = @(X)getIfdl(meta,X,0);
    data.tfDIfdl = @(r,X)getIfdl(meta,X,r);
    data.tfPhise = @(X)getPhise(meta,X,0);
    data.tfDPhise = @(r,X)getPhise(meta,X,r);
    data.tfThetass = @(X)getThetass(meta,X,0);
    data.tfDThetass = @(r,X)getThetass(meta,X,r);
    data.tfEta = @(X)getEta(meta,X,0);
    data.tfDEta = @(r,X)getEta(meta,X,r);
    data.tfPhie = @(X)getPhie(meta,X);
    data.tfVcell = @()getVcell(meta);
    data.meta = meta;

    function C = solveTwoInertLayers
        % Allocate storage for linear system of the form Ac = b.
        % Some entries of the A and b matricies are constant with frequency.
        A = zeros(8); 
        A(3,1) = +1;
        A(3,4) = -1;
        A(5,3) = +1;
        A(5,5) = -1;
        A(5,7) = -1;
        b = zeros(8,1);
        b(1) = -1/p.kappad/p.psi/p.T;
        
        % Allocate storage for coefficients.
        C = zeros(ns,8);
        
        for i = 1:ns
            % Update A matrix.
            LS  = LambdaS(i);  expnLS = expnLambdaS(i);
            LD  = LambdaD(i);  expnLD = expnLambdaD(i);
            L1P = Lambda1P(i); expnL1P = expnLambda1P(i);
            L2P = Lambda2P(i); expnL2P = expnLambda2P(i);
            l1p = lambda1p(i);
            l2p = lambda2p(i);
            A(1,1) = LD*expnLD;
            A(1,2) = -LD;
            A(2,1) = p.kappad*LD;
            A(2,2) = -p.kappad*LD*expnLD;
            A(2,3) = -p.kappas*LS*expnLS;
            A(2,4) = p.kappas*LS;
            A(3,2) = expnLD;
            A(3,3) = -expnLS;
            A(4,3) = p.kappas*LS;
            A(4,4) = -p.kappas*LS*expnLS;
            A(4,5) = p.kappap*L1P;
            A(4,6) = -p.kappap*L1P*expnL1P;
            A(4,7) = p.kappap*L2P;
            A(4,8) = -p.kappap*L2P*expnL2P;
            A(5,4) = expnLS;
            A(5,6) = -expnL1P;
            A(5,8) = -expnL2P;
            A(6,5) = l1p;
            A(6,6) = -l1p*expnL1P;
            A(6,7) = l2p;
            A(6,8) = -l2p*expnL2P;
            A(7,5) = L1P*expnL1P;
            A(7,6) = -L1P;
            A(7,7) = L2P*expnL2P;
            A(7,8) = -L2P;
            A(8,5) = expnL1P*(L1P^3);
            A(8,6) = -(L1P^3);
            A(8,7) = expnL2P*(L2P^3);
            A(8,8) = -(L2P^3);
        
            % Update b matrix.
            b(6) = mu1p(i)/p.kappap;
            b(8) = -mu1p(i)/p.sigmap;
        
            % Solve linear system for coefficients.
            c = A\b;
        
            % Store result.
            C(i,:) = c(:).';
        end % i=1..n
    end
    
    function C = solveSingleInertLayer
        % Allocate storage for linear system of the form Ac = b.
        % Some entries of the A and b matricies are constant with frequency.
        A = zeros(6); 
        A(2,1) = +1;
        A(2,3) = -1;
        A(2,5) = -1;
        b = zeros(6,1);
        b(1) = -1/p.kappas/p.psi/p.T;
        
        % Allocate storage for coefficients.
        C = zeros(ns,6);
        
        for i = 1:ns
            % Update A matrix.
            LS  = LambdaS(i);  expnLS = expnLambdaS(i);
            L1P = Lambda1P(i); expnL1P = expnLambda1P(i);
            L2P = Lambda2P(i); expnL2P = expnLambda2P(i);
            l1p = lambda1p(i);
            l2p = lambda2p(i);
            A(1,1) = LS*expnLS;
            A(1,2) = -LS;
            A(2,2) = expnLS;
            A(2,4) = -expnL1P;
            A(2,6) = -expnL2P;
            A(3,1) = LS*p.kappas;
            A(3,2) = -LS*p.kappas*expnLS;
            A(3,3) = L1P*p.kappap;
            A(3,4) = -L1P*p.kappap*expnL1P;
            A(3,5) = L2P*p.kappap;
            A(3,6) = -L2P*p.kappap*expnL2P;
            A(4,3) = l1p;
            A(4,4) = -l1p*expnL1P;
            A(4,5) = l2p;
            A(4,6) = -l2p*expnL2P;
            A(5,3) = L1P*expnL1P;
            A(5,4) = -L1P;
            A(5,5) = L2P*expnL2P;
            A(5,6) = -L2P;
            A(6,3) = expnL1P*(L1P^3);
            A(6,4) = -(L1P^3);
            A(6,5) = expnL2P*(L2P^3);
            A(6,6) = -(L2P^3);
        
            % Update b matrix.
            b(4) = mu1p(i)/p.kappap;
            b(6) = -mu1p(i)/p.sigmap;
        
            % Solve linear system for coefficients.
            c = A\b;
        
            % Store result.
            C(i,:) = c(:).';
        end % i=1..n
    end
end

function Thetae = getThetae(meta,X,r)
    [Z, bins] = meta.param.split(X);
    Thetae = zeros(meta.ns,length(X));

    if any(bins.n)
        if isfield(meta,'c1d')
            LD = meta.LambdaD;
            LDr = LD.^r;
            negLDr = LDr.*(-1).^r;
            Thetae(:,bins.n) = ...
                meta.c1d.*LDr.*exp(LD.*(Z.n-1)) + meta.c2d.*negLDr.*exp(-LD.*Z.n);
        else
            LS = meta.LambdaS;
            LSr = LS.^r;
            negLSr = LSr.*(-1).^r;
            Thetae(:,bins.n) = ...
                meta.c1s.*LSr.*exp(LS.*(Z.n-1)) + meta.c2s.*negLSr.*exp(-LS.*Z.n);
        end
    end

    if any(bins.d)
        LD = meta.LambdaD;
        LDr = LD.^r;
        negLDr = LDr.*(-1).^r;
        Thetae(:,bins.d) = ...
            meta.c1d.*LDr.*exp(LD.*(Z.d-1)) + meta.c2d.*negLDr.*exp(-LD.*Z.d);
    end
    
    if any(bins.s)
        LS = meta.LambdaS;
        LSr = LS.^r;
        negLSr = LSr.*(-1).^r;
        Thetae(:,bins.s) = ...
            meta.c1s.*LSr.*exp(LS.*(Z.s-1)) + meta.c2s.*negLSr.*exp(-LS.*Z.s);
    end

    if any(bins.p)
        L1P = meta.Lambda1P;
        L1Pr = L1P.^r;
        negL1Pr = L1Pr.*(-1).^r;
        L2P = meta.Lambda2P;
        L2Pr = L2P.^r;
        negL2Pr = L2Pr.*(-1).^r;
        Thetae(:,bins.p) = ...
            meta.c1p.*L1Pr.*exp(L1P.*(Z.p-1)) + meta.c2p.*negL1Pr.*exp(-L1P.*Z.p) + ...
            meta.c3p.*L2Pr.*exp(L2P.*(Z.p-1)) + meta.c4p.*negL2Pr.*exp(-L2P.*Z.p);
        Thetae(:,bins.p) = Thetae(:,bins.p).*(-1).^r; % dz/dx = -1
    end
end

function Ifdl = getIfdl(meta,X,r)
    [Z, bins] = meta.param.split(X);
    Ifdl = zeros(meta.ns,length(X));

    if any(bins.n)
        if r == 0
            Ifdl(:,bins.n) = 1;  % Ifdl(n)=Iapp!
        else
            Ifdl(:,bins.n) = NaN;
        end
    end

    if any(bins.d)
        Ifdl(:,bins.d) = 0;  % Ifdl(d)=0!
    end
    
    if any(bins.s)
        Ifdl(:,bins.s) = 0;  % Ifdl(s)=0!
    end

    if any(bins.p)
        L1P = meta.Lambda1P;
        L1Pr = L1P.^r;
        negL1Pr = L1Pr.*(-1).^r;
        L2P = meta.Lambda2P;
        L2Pr = L2P.^r;
        negL2Pr = L2Pr.*(-1).^r;
        Ifdl(:,bins.p) = meta.j1p.*L1Pr.*exp(L1P.*(Z.p-1)) ...
            + meta.j2p.*negL1Pr.*exp(-L1P.*Z.p) ...
            + meta.j3p.*L2Pr.*exp(L2P.*(Z.p-1)) ...
            + meta.j4p.*negL2Pr.*exp(-L2P.*Z.p);
        Ifdl(:,bins.p) = Ifdl(:,bins.p).*(-1).^r; % dz/dx = -1
    end
end

function Phise = getPhise(meta,X,r)
    [~, bins] = meta.param.split(X);
    Phise = zeros(meta.ns,length(X));

    if any(bins.n)
        if r == 0
            Phise(:,bins.n) = meta.param.Zsen;
        else
            Phise(:,bins.n) = NaN;
        end
    end

    if any(bins.d)
        Phise(:,bins.d) = NaN;
    end

    if any(bins.s)
        Phise(:,bins.s) = NaN;
    end

    if any(bins.p)
        Phise(:,bins.p) = meta.param.Zsep.*getIfdl(meta,X(bins.p),r);
    end
end

function Thetass = getThetass(meta,X,r)
    [~, bins] = meta.param.split(X);
    Thetass = zeros(meta.ns,length(X));

    if any(bins.n|bins.d|bins.s)
        Thetass(:,bins.n|bins.d|bins.s) = NaN;
    end

    if any(bins.p)
        Thetass(:,bins.p) = ...
            (meta.param.Zsp.*meta.divFactorIfp./meta.param.dUocpp)...
            .*getIfdl(meta,X(bins.p),r);
    end
end

function Eta = getEta(meta,X,r)
    [~, bins] = meta.param.split(X);
    Eta = zeros(meta.ns,length(X));

    if any(bins.n)
        if r == 0
            Rctn = meta.param.Rctn;
            divFactorIfn = meta.divFactorIfn;
            Eta(:,bins.n) = divFactorIfn.*Rctn.*ones(1,sum(bins.n));
        else
            Eta(:,bins.n) = NaN;
        end
    end

    if any(bins.d|bins.s)
        Eta(:,bins.d|bins.s) = NaN;
    end

    if any(bins.p)
        Rctp = meta.param.Rctp;
        divFactorIfp = meta.divFactorIfp;
        Eta(:,bins.p) = ...
            (Rctp.*divFactorIfp).*getIfdl(meta,X(bins.p),r);
    end
end

function [Phie, data] = getPhie(meta,X)
    [Z, bins] = meta.param.split(X);
    Phie1 = zeros(meta.ns,length(X));

    kappad = meta.param.kappad;
    kappas = meta.param.kappas;
    kappap = meta.param.kappap;

    if any(bins.n)
        Phie1(:,bins.n) = 0;
    end

    if any(bins.d)
        Phie1(:,bins.d) = -( Z.d/kappad );
    end

    if any(bins.s)
        Phie1(:,bins.s) = -( 1/kappad + Z.s/kappas );
    end

    if any(bins.p)
        L1P = meta.Lambda1P;
        L1Psq = L1P.^2;
        L2P = meta.Lambda2P;
        L2Psq = L2P.^2;
        Phie1(:,bins.p) = -1/kappad - 1/kappas - (1-Z.p)/kappap ...
            - (meta.j1p./L1Psq./kappap).*(exp(L1P.*(Z.p-1))-L1P.*(Z.p-1)-1) ...
            - (meta.j2p./L1Psq./kappap).*(exp(-L1P.*Z.p)+L1P.*(Z.p-1).*exp(-L1P)-exp(-L1P)) ...
            - (meta.j3p./L2Psq./kappap).*(exp(L2P.*(Z.p-1))-L2P.*(Z.p-1)-1) ...
            - (meta.j4p./L2Psq./kappap).*(exp(-L2P.*Z.p)+L2P.*(Z.p-1).*exp(-L2P)-exp(-L2P));
    end

    Thetae = getThetae(meta,X,0);
    Thetae0 = getThetae(meta,0,0);
    Phie2 = -meta.param.kD*meta.T*(Thetae-Thetae0);
    Phie = Phie1 + Phie2;
    data.Phie1 = Phie1;
    data.Phie2 = Phie2;
end

function [Vcell, data] = getVcell(meta)
    Phise = getPhise(meta,[meta.param.xlim.min meta.param.xlim.max],0);
    Phise0 = Phise(:,1);
    PhiseL = Phise(:,2);
    PhieL = getPhie(meta,meta.param.xlim.max);
    Vcell = -(PhiseL + PhieL - Phise0);
    data.Zneg = Phise0;
    data.Zpos = -PhiseL;
    data.Zel = -PhieL;
end


% Second-Harmonic ---------------------------------------------------------

function data = tfLMB22(lindata,layerMode)
    p = lindata.meta.param;
    S = lindata.meta.S;
    ns = length(S);

    mu22p = p.kD/p.psi/p.kappap./p.Zse2p - ((3600*p.qep/p.psi/p.T/p.kappap).*2*S);
    tau12p = (1/p.sigmap + 1/p.kappap)./p.Zse2p - mu22p; 
    tau22p = ((3600*p.qep/p.psi/p.T/p.kappap).*2*S).*(1/p.sigmap + 1/p.kappap)./p.Zse2p;
    LambdaS2 = (3600*p.qes*2*S/p.kappas/p.psi/p.T);
    LambdaD2 = (3600*p.qed*2*S/p.kappad/p.psi/p.T);

    % BVP solver options.
    opts = bvpset('AbsTol',1e-6,'RelTol',1e-3);  % defaults

    % Mesh points for BVP solver.
    zmesh = linspace(0,1,10);

    % Points at which to evalulate BVP solution for Thetae22 
    % (will be used for numerical integration / interpolation later).
    zsoln = linspace(0,1,1000);

    % Run BVP solver (see function definitions below!)
    if strcmpi(layerMode,'off')
        soln = solveTwoIntertLayers;
    else
        soln = solveSingleInertLayer;
    end

    % Build metadata structure.
    meta = getNL22(lindata,soln.xsoln);
    meta.divFactorIfp = p.Zdl2p./(p.Rctp+p.Zs2p+p.Zdl2p);
    meta.divFactorIdlp = (p.Rctp+p.Zs2p)./(p.Rctp+p.Zs2p+p.Zdl2p);
    meta.divFactorIfn = p.Zdl2n./(p.Rctn+p.Zdl2n);
    meta.divFactorIdln = p.Rctn./(p.Rctn+p.Zdl2n);
    meta.Thetae22Init = soln.Thetae22Init;
    meta.Thetae22 = soln.Thetae22;
    meta.d1Thetae22 = soln.d1Thetae22;
    meta.d2Thetae22 = soln.d2Thetae22;
    meta.i2Thetae22 = soln.i2Thetae22;
    meta.xsoln = soln.xsoln;
    meta.indneg = find(soln.xsoln==p.xlim.min);
    meta.indposlwr = find(soln.xsoln==p.xlim.sp);
    meta.indposupr = find(soln.xsoln==p.xlim.max);
    meta.S = S;
    meta.ns = length(S);
    meta.TdegC = p.TdegC;
    meta.T = p.T;
    meta.R = p.R;
    meta.F = p.F;
    meta.param = p;
    meta.aux.mu22p = mu22p;
    meta.aux.tau12p = tau12p;
    meta.aux.tau22p = tau22p;
    meta.aux.LambdaS2 = LambdaS2;
    meta.aux.LambdaD2 = LambdaD2;

    % Build output data structure.
    data.tfThetaeInit = @(X)getThetae22Init(meta,X);
    data.tfThetae = @(X)getThetae22(meta,X);
    data.tfIfdl = @(X)getIfdl22(meta,X,true);
    data.tfPhise = @(X)getPhise22(meta,X,true);
    data.tfPhie = @(X)getPhie22(meta,X,true);
    data.tfThetass = @(X)getThetass22(meta,X,true);
    data.tfEta = @(X)getEta22(meta,X,true);
    data.tfVcell = @()getVcell22(meta);
    data.meta = meta;

    function soln = solveSingleInertLayer
        ind.Thetaep = 1;
        ind.d1Thetaep = 2;
        ind.d2Thetaep = 3;
        ind.d3Thetaep = 4;
        ind.Thetaes = 5;
        ind.d1Thetaes = 6;

        % In the BVP function, we need to evalulate d2g/dx2 at the present
        % x-location in the positive electrode. We evalulate d2g/dx2 at the
        % points given below and linearly interpolate inside of the BVP
        % function to make evalulation fast.
        xInterpPos = linspace(1,2,1000);
    
        % Collect quantities needed to solve BVP
        nl22 = getNL22(lindata,[1 2]);
        dg1 = nl22.dg(:,1);  % dg/dx at x=1
        dg2 = nl22.dg(:,2);  % dg/dx at x=2
        nl22 = getNL22(lindata,xInterpPos);
        d2gInterp = nl22.d2g;
    
        % Linear system of the form Ac=b for generating initial guess for 
        % Thetae22. c is a vector of polynomial coefficients. See
        % ThetaeInitial() below for more details.
        A = zeros(5);
        A(1,:) = [1 -1 -1 -1 -1];
        A(2,:) = [2*p.kappas [-4 -4 -2 -1]*p.kappap];
        A(4,2:end) = [32 12 4 1];
        A(5,2:3) = [48 6];
        b = zeros(5,1);
        
        % Solve 4th-order BVP in Thetae at each frequency.
        % Allocate storage for solution.
        Thetae22sInit = zeros(ns,length(zsoln));
        Thetae22pInit = zeros(ns,length(zsoln));
        Thetae22s = zeros(ns,length(zsoln));
        Thetae22p = zeros(ns,length(zsoln));
        d2Thetae22p = zeros(ns,length(zsoln));  % only available in pos
        for k = 1:ns
            tau1p = tau12p(k);
            tau2p = tau22p(k);
            mu2p = mu22p(k);
    
            % Update linear system.
            A(3,2:end) = [24+4*mu2p 6+3*mu2p 2*mu2p mu2p];
            b(3) = dg1(k);
            b(5) = dg2(k);
    
            % Solve for polynomial coefficients.
            C = A\b;
    
            % Initialize BVP solver.
            sol = bvpinit(zmesh,@(z)ThetaeInitial(C,z));
            Thetae22pInit(k,:) = interp1(sol.x,sol.y(ind.Thetaep,:),zsoln,'linear','extrap');
            Thetae22sInit(k,:) = interp1(sol.x,sol.y(ind.Thetaes,:),zsoln,'linear','extrap');
    
            % Solve BVP.
            sol = bvp5c(@ThetaeODE, @ThetaeBC, sol, opts);
            Thetae = deval(sol,zsoln,[ind.Thetaep ind.d2Thetaep ind.Thetaes]);
            Thetae22p(k,:) = Thetae(1,:);
            d2Thetae22p(k,:) = Thetae(2,:);
            Thetae22s(k,:) = Thetae(3,:);
        end
    
        soln.xsoln = [zsoln(1:end-1) 1+zsoln]; % x=1+z in pos
        soln.Thetae22Init = [Thetae22sInit(:,1:end-1) Thetae22pInit];
        soln.Thetae22 = [Thetae22s(:,1:end-1) Thetae22p];
        soln.d2Thetae22 = [nan(ns,length(zsoln)-1) d2Thetae22p];
        soln.iThetae22p = cumtrapz(zsoln,Thetae22p.').';
        soln.i2Thetae22p = cumtrapz(zsoln,iThetae22p.').';
        soln.i2Thetae22 = [nan(ns,length(zsoln)-1) i2Thetae22p];

        function dydx = ThetaeODE(z,y)
            % Compute derivative of the y-vector.
            % Let: y = [
            %         Thetae(p)
            %         Thetae(p)'
            %         Thetae(p)''
            %         Thetae(p)'''
            %         ---------
            %         Thetae(s)
            %         Thetae(s)'
            %      ]
            % where prime denotes the partial w/r/t x(tilde).
            % Then: y' = [
            %         y(2)
            %         y(3)
            %         y(4)
            %         tau1*y(3) - tau2*y(1) + g(x(tilde))
            %         ---------
            %         y(6)
            %         [(3600*qe*2*s)/(kappa*psi*T)]*y(5)
            %       ]
    
            d2gx = fastinterp(xInterpPos,d2gInterp(k,:),1+z);  % x=1+z in pos
            dydx = [
                % --- pos ---
                y(ind.d1Thetaep)
                y(ind.d2Thetaep)
                y(ind.d3Thetaep)
                tau1p*y(ind.d2Thetaep) - tau2p*y(ind.Thetaep) + d2gx
                % --- sep ---
                y(ind.d1Thetaes)
                LambdaS2(k)*y(ind.Thetaes)
            ];
        end
    
        function residual = ThetaeBC(ya,yb)
            % Compute boundary-condition residuals.
            % See ThetaeODE() for the components of the y-vector.
    
            residual = [
                % --- sep ---
                ya(ind.d1Thetaes)                 % 1. zero salt flux at x=0
                % --- interface ---
                ya(ind.Thetaep)-yb(ind.Thetaes)   % 2. continuity of salt conc. at x=1
                p.kappap*ya(ind.d1Thetaep)-p.kappas*yb(ind.d1Thetaes) % 3. continuity of salt flux at x=1
                % --- pos ---
                ya(ind.d3Thetaep)+mu2p*ya(ind.d1Thetaep)-dg1(k) % 4. zero solid/salt current at x=1+
                yb(ind.d1Thetaep)                 % 5. zero salt flux at x=2
                yb(ind.d3Thetaep)-dg2(k)          % 6. zero solid/salt current at x=2
            ];
        end
    
        function yinit = ThetaeInitial(C,Z)
            % Generate initial guess for the Thetae y-vector. 
            % Polynomial approximation.
            % In the separator, we assume:
            %    Thetae = a*X^2 + b*X + c
            % In the porous electrode, we assume:
            %    Thetae = d*X^4 + e*X^3 + f*X^2 + h*X + p
            % When applying the boundary conditions, we find that b=0. We also
            % assume c=0 and p=0. C is a vector of the other polynomial
            % coefficients chosen to satisfy the remaining five boundary
            % conditions:
            %    C = [a d e f h]
            % Specifically, C is the solution to the following linear system:
            %   
            %    1        -1        -1        -1        -1                 0
            %    2ks      -4kp      -3kp      -2kp      -kp                0
            %    0        24+4mu2p  6+2mu2p   2mu2p     mu2p    *   C  =   dg1
            %    0        32        12        4         1                  0
            %    0        48        6         0         0                  dg2
            %
            % where ks=kappas, kp=kappap, dg1 = dg/dx at x=1, and dg2 = dg/dx 
            % at x=2.
        
            % Unpack polynomial coefficients.
            ca = C(1);
            cb = 0;
            cc = 0;
            cd = C(2);
            ce = C(3);
            cf = C(4);
            ch = C(5);
            cp = 0;
    
            Xp = 1+Z;
            yinit = [
                % --- pos ---
                cd*Xp.^4 + ce*Xp.^3 + cf*Xp.^2 + ch*Xp + cp
                4*cd*Xp.^3 + 3*ce*Xp.^2 + 2*cf*Xp + ch
                12*cd*Xp.^2 + 6*ce*Xp + 2*cf
                24*cd*Xp + 6*ce
                % --- sep ---
                ca*Z.^2 + cb*Z + cc
                2*ca*Z + cb
            ];
        end
    end % solveSingleInertLayer()

    function soln = solveTwoIntertLayers
        ind.Thetaep = 1;
        ind.d1Thetaep = 2;
        ind.d2Thetaep = 3;
        ind.d3Thetaep = 4;
        ind.Thetaes = 5;
        ind.d1Thetaes = 6;
        ind.Thetaed = 7;
        ind.d1Thetaed = 8;

        % In the BVP function, we need to evalulate d2g/dx2 at the present
        % x-location in the positive electrode. We evalulate d2g/dx2 at the
        % points given below and linearly interpolate inside of the BVP
        % function to make evalulation fast.
        xInterpPos = linspace(2,3,1000);
    
        % Collect quantities needed to solve BVP
        nl22 = getNL22(lindata,[2 3]);
        dg1 = nl22.dg(:,1);  % dg/dx at x=2
        dg2 = nl22.dg(:,2);  % dg/dx at x=3
        nl22 = getNL22(lindata,xInterpPos);
        d2gInterp = nl22.d2g;
        
        % Solve 4th-order BVP in Thetae at each frequency.
        % Allocate storage for solution.
        Thetae22pInit = zeros(ns,length(zsoln));
        Thetae22sInit = zeros(ns,length(zsoln));
        Thetae22dInit = zeros(ns,length(zsoln));
        Thetae22p = zeros(ns,length(zsoln));
        Thetae22s = zeros(ns,length(zsoln));
        Thetae22d = zeros(ns,length(zsoln));
        d1Thetae22p = zeros(ns,length(zsoln));
        d1Thetae22s = zeros(ns,length(zsoln));
        d1Thetae22d = zeros(ns,length(zsoln));
        d2Thetae22p = zeros(ns,length(zsoln));  % only available in pos
        for k = 1:ns
            tau1p = tau12p(k);
            tau2p = tau22p(k);
            mu2p = mu22p(k);
    
            % Initialize BVP solver.
            sol = bvpinit(zmesh,@(z)ThetaeInitial([],z));
            Thetae22pInit(k,:) = interp1(sol.x,sol.y(ind.Thetaep,:),zsoln,'linear','extrap');
            Thetae22sInit(k,:) = interp1(sol.x,sol.y(ind.Thetaes,:),zsoln,'linear','extrap');
            Thetae22dInit(k,:) = interp1(sol.x,sol.y(ind.Thetaed,:),zsoln,'linear','extrap');
    
            % Solve BVP.
            sol = bvp5c(@ThetaeODE, @ThetaeBC, sol, opts);
            Thetae = deval(sol,zsoln,[
                ind.Thetaep
                ind.d1Thetaep
                ind.d2Thetaep
                ind.Thetaes
                ind.d1Thetaes
                ind.Thetaed
                ind.d1Thetaed
            ].');
            Thetae22p(k,:) = Thetae(1,:);
            d1Thetae22p(k,:) = Thetae(2,:);
            d2Thetae22p(k,:) = Thetae(3,:);
            Thetae22s(k,:) = Thetae(4,:);
            d1Thetae22s(k,:) = Thetae(5,:);
            Thetae22d(k,:) = Thetae(6,:);
            d1Thetae22d(k,:) = Thetae(7,:);
        end
    
        soln.xsoln = [zsoln(1:end-1) 1+zsoln(1:end-1) 2+zsoln]; % x=1+z in sep; x=2+z in pos
        soln.Thetae22Init = [Thetae22dInit(:,1:end-1) Thetae22sInit(:,1:end-1) Thetae22pInit];
        soln.Thetae22 = [Thetae22d(:,1:end-1) Thetae22s(:,1:end-1) Thetae22p];
        soln.d1Thetae22 = [d1Thetae22d(:,1:end-1) d1Thetae22s(:,1:end-1) d1Thetae22p];
        soln.d2Thetae22 = [nan(ns,2*(length(zsoln)-1)) d2Thetae22p];
        soln.i1Thetae22p = cumtrapz(zsoln,Thetae22p.').';
        soln.i2Thetae22p = cumtrapz(zsoln,soln.i1Thetae22p.').';
        soln.i2Thetae22 = [nan(ns,2*(length(zsoln)-1)) soln.i2Thetae22p];

        function dydx = ThetaeODE(z,y)
            % Compute derivative of the y-vector.
            % Let: y = [
            %         Thetae(p)
            %         Thetae(p)'
            %         Thetae(p)''
            %         Thetae(p)'''
            %         ---------
            %         Thetae(s)
            %         Thetae(s)'
            %         ---------
            %         Thetae(d)
            %         Thetae(d)'
            %      ]
            % where prime denotes the partial w/r/t x(tilde).
            % Then: y' = [
            %         y(2)
            %         y(3)
            %         y(4)
            %         tau1*y(3) - tau2*y(1) + g(x(tilde))
            %         ---------
            %         y(6)
            %         [(3600*qes*2*s)/(kappas*psi*T)]*y(5)
            %         ---------
            %         y(8)
            %         [(3600*qed*2*s)/(kappad*psi*T)]*y(7)
            %       ]
    
            d2gx = fastinterp(xInterpPos,d2gInterp(k,:),2+z); % x=2+z in pos
            dydx = [
                % --- pos ---
                y(ind.d1Thetaep)
                y(ind.d2Thetaep)
                y(ind.d3Thetaep)
                tau1p*y(ind.d2Thetaep) - tau2p*y(ind.Thetaep) + d2gx
                % --- sep ---
                y(ind.d1Thetaes)
                LambdaS2(k)*y(ind.Thetaes)
                % --- dll ---
                y(ind.d1Thetaed)
                LambdaD2(k)*y(ind.Thetaed)
            ];
        end
    
        function residual = ThetaeBC(ya,yb)
            % Compute boundary-condition residuals.
            % See ThetaeODE() for the components of the y-vector.
    
            residual = [
                % --- dll ---
                ya(ind.d1Thetaed)                 % 1. zero salt flux at x=0
                % --- dll-sep interface ---
                ya(ind.Thetaes)-yb(ind.Thetaed)   % 2. continuity of salt conc. at x=1
                p.kappas*ya(ind.d1Thetaes)-p.kappad*yb(ind.d1Thetaed) % 3. continuity of salt flux at x=1
                % --- sep-pos interface ---
                ya(ind.Thetaep)-yb(ind.Thetaes)   % 4. continuity of salt conc. at x=2
                p.kappap*ya(ind.d1Thetaep)-p.kappas*yb(ind.d1Thetaes) % 5. continuity of salt flux at x=2
                % --- pos ---
                ya(ind.d3Thetaep)+mu2p*ya(ind.d1Thetaep)-dg1(k) % 6. zero solid/salt current at x=2+
                yb(ind.d1Thetaep)                 % 7. zero salt flux at x=3
                yb(ind.d3Thetaep)-dg2(k)          % 8. zero solid/salt current at x=3
            ];
        end
    
        function yinit = ThetaeInitial(~,Z)
            % Generate initial guess for the Thetae y-vector.

            % Trivial solution. (Does not satisfy some of the pos BCs.)
            yinit = zeros(length(fieldnames(ind)),length(Z));
        end
    end % solveTwoInertLayers()
end

function Thetae = getThetae22Init(meta,X)
    Thetae = interp1(meta.xsoln,meta.Thetae22Init.',X).';
end

function Thetae = getThetae22(meta,X)
    Thetae = interp1(meta.xsoln,meta.Thetae22.',X).';
end

function Thetae = getd1Thetae22(meta,X)
    Thetae = interp1(meta.xsoln,meta.d1Thetae22.',X).';
end

function Ifdl = getIfdl22(meta,X,delinearize)
    [~, bins] = meta.param.split(X);
    Ifdl = zeros(meta.ns,length(X));

    if any(bins.n)
        if delinearize
            Ifdl(:,bins.n) = 0;  % ifdl(n)=0!
        else
            Ifdl(:,bins.n) = -meta.DeltaIfdl(:,meta.indneg);  % Ifdl(n)=-DeltaIdl11(n)!
        end
    end
    
    if any(bins.s)
        Ifdl(:,bins.s) = 0;  % Ifdl(s)=0!
    end

    if any(bins.p)
        Thetae = interp1(meta.xsoln,meta.Thetae22.',X(bins.p)).';
        d2Thetae = interp1(meta.xsoln,meta.d2Thetae22.',X(bins.p)).';
        qe = meta.param.qep;
        kappa = meta.param.kappap;
        psi = meta.param.psi;
        T = meta.param.T;
        if delinearize
            Ifdl(:,bins.p) = 3600*qe*2*meta.S.*Thetae - kappa*psi*T*d2Thetae;
        else
            DeltaIfdl11 = interp1(meta.xsoln,meta.DeltaIfdl.',X(bins.p)).';
            Ifdl(:,bins.p) = 3600*qe*2*meta.S.*Thetae ...
               - kappa*psi*T*d2Thetae - DeltaIfdl11;
        end
    end
end

function If = getIf22(meta,X,delinearize)
    [~, bins] = meta.param.split(X);
    If = zeros(meta.ns,length(X));

    if any(bins.n)
        if delinearize
            If(:,bins.n) = 0;  % if(n)=0!
        else
            If(:,bins.n) = -meta.divFactorIfn*meta.DeltaIfdl(:,meta.indneg);
        end
    end
    
    if any(bins.s)
        If(:,bins.s) = 0;  % If(s)=0!
    end

    if any(bins.p)
        Ifdl = getIfdl22(meta,X(bins.p),false);
        If = meta.divFactorIfp.*Ifdl;
        if delinearize
            If = If + interp1(meta.xsoln,meta.DeltaIf.',X(bins.p)).';
        end
        If(:,bins.p) = If;
    end
end

function Thetass = getThetass22(meta,X,delinearize)
    [~, bins] = meta.param.split(X);
    Thetass = zeros(meta.ns,length(X));

    if any(bins.s|bins.n)
        Thetass(:,bins.s|bins.n) = NaN;
    end

    if any(bins.p)
        If = getIf22(meta,X(bins.p),delinearize);
        Thetass(:,bins.p) = If.*meta.param.Zs2p/meta.param.dUocpp;
    end
end

function Phise = getPhise22(meta,X,delinearize)
    [~, bins] = meta.param.split(X);
    Phise = zeros(meta.ns,length(X));

    if any(bins.n)
        Zse2 = meta.param.Zse2n;
        Phise(:,bins.n) = -Zse2.*meta.DeltaIfdl(:,meta.indneg);
        if delinearize
            Phise(:,bins.n) = Phise(:,bins.n) ... 
                + meta.DeltaPhis(:,meta.indneg) ...
                - meta.DeltaPhie(:,meta.indneg);
        end
    end

    if any(bins.s)
        Phise(:,bins.s) = NaN;
    end

    if any(bins.p)
        Ifdl22 = getIfdl22(meta,X(bins.p),false);
        DeltaIf11 = interp1(meta.xsoln,meta.DeltaIf.',X(bins.p)).';
        Zse2p = meta.param.Zse2p;
        Zs2p = meta.param.Zs2p;
        divFactorIfp = meta.divFactorIfp;
        Phise(:,bins.p) = Zse2p.*Ifdl22+divFactorIfp.*Zs2p.*DeltaIf11;
        if delinearize
            DeltaPhise11 = interp1(meta.xsoln,meta.DeltaPhise.',X(bins.p)).';
            Phise(:,bins.p) = Phise(:,bins.p) + DeltaPhise11;
        end
    end
end

function Eta = getEta22(meta,X,delinearize)
    [~, bins] = meta.param.split(X);
    Eta = zeros(meta.ns,length(X));

    if any(bins.n)
        Zse2 = meta.param.Zse2n;
        Rf = meta.param.Rfn;
        Eta(:,bins.n) = -(Zse2+Rf).*meta.DeltaIfdl(:,meta.indneg);
        if delinearize
            Eta(:,bins.n) = Eta(:,bins.n) + meta.DeltaEta(:,meta.indneg);
        end
    end

    if any(bins.s)
        Eta(:,bins.s) = NaN;
    end

    if any(bins.p)
        Rf = meta.param.Rfp;
        Ifdl22 = getIfdl22(meta,X(bins.p),false);
        DeltaIf11 = interp1(meta.xsoln,meta.DeltaIf.',X(bins.p)).';
        Zse2p = meta.param.Zse2p;
        Zs2p = meta.param.Zs2p;
        divFactorIfp = meta.divFactorIfp;
        if22 = meta.divFactorIfp.*Ifdl22 ...
            + interp1(meta.xsoln,meta.DeltaIf.',X(bins.p)).';
        Eta(:,bins.p) = Zse2p.*Ifdl22+divFactorIfp.*Zs2p.*DeltaIf11 ...
            - Ifdl22*Rf - meta.param.Zs2p.*if22;
        if delinearize
            DeltaEta11 = interp1(meta.xsoln,meta.DeltaEta.',X(bins.p)).';
            Eta(:,bins.p) = Eta(:,bins.p) + DeltaEta11;
        end
    end
end

function [Phie, data] = getPhie22(meta,X,delinearize)
    [~, bins] = meta.param.split(X);
    Phie1 = zeros(meta.ns,length(X));
    
    Thetae = getThetae22(meta,X);
    Thetae0 = getThetae22(meta,0);
    Phie2 = -meta.param.kD*meta.T*(Thetae-Thetae0);

    if any(bins.n|bins.d|bins.s)
        Phie1(:,bins.n|bins.d|bins.s) = 0;
    end

    if any(bins.p)
        kappap = meta.param.kappap;
        qep = meta.param.qep;
        psi = meta.param.psi;
        T = meta.T;
        S = meta.S;
        ThetaeSP = getThetae22(meta,meta.param.xlim.sp);
        d1ThetaeSP = getd1Thetae22(meta,meta.param.xlim.sp);
        i2Thetae22 = interp1(meta.xsoln,meta.i2Thetae22.',X(bins.p)).';
        Phie1(:,bins.p) = psi*T*(Thetae(:,bins.p)-ThetaeSP) ...
            - (3600*qep*2*S/kappap).*i2Thetae22 ...
            - psi*T*d1ThetaeSP*(X(bins.p)-meta.param.xlim.sp);
    end

    Phie = Phie1 + Phie2;
    if delinearize
        Phie = Phie + interp1(meta.xsoln,meta.DeltaPhie.',X).';
    end

    data.Phie1 = Phie1;
    data.Phie2 = Phie2;
end

function [Vcell, data] = getVcell22(meta)
    Phise = getPhise22(meta,[meta.param.xlim.min meta.param.xlim.max],true);
    Phise0 = Phise(:,1);
    PhiseL = Phise(:,2);
    PhieL = getPhie22(meta,meta.param.xlim.max,true);
    Vcell = -(PhiseL + PhieL - Phise0);
    data.Zneg = Phise0;
    data.Zpos = -PhiseL;
    data.Zel = -PhieL;
end


% Utilities ---------------------------------------------------------------

function y0 = fastinterp(x,y,x0)
    %FASTINTERP Fast linear interpolation. Single-point.

    k = find(x>=x0,1,'first');
    if x0==x(k)
        y0 = y(k);
        return;
    end
    x1 = x(k-1);
    x2 = x(k);
    y1 = y(k-1);
    y2 = y(k);
    m = (y2-y1)/(x2-x1);
    y0 = y1 + m*(x0-x1);
end

function data = getNL22(lindata,X)
    %GETNL22 Compute nonlinear correction terms for second harmonic.

    meta = lindata.meta;
    [~, bins] = meta.param.split(X);

    psi = meta.param.psi;
    kD = meta.param.kD;
    T = meta.T;
    R = meta.R;
    F = meta.F;
    kappap = meta.param.kappap;
    Rfp = meta.param.Rfp;
    Rct1p = meta.param.Rctp;
    Rct2invp = meta.param.Rct2invp;
    d2Uocpp = meta.param.d2Uocpp;
    Rfn = meta.param.Rfn;
    Rct1n = meta.param.Rctn;
    Rct2invn = meta.param.Rct2invn;
    Zdl2p = meta.param.Zdl2p;
    Zs2p = meta.param.Zs2p;
    Zse2p = meta.param.Zse2p;
    Zdl2n = meta.param.Zdl2n;

    % Note: computing g(x).
    % We need thetass11^2, thetae11^2, eta11^2 and first/second
    % derivatives thereof w/r/t x(tilde) to compute Delta(If11), 
    % Delta(If+dl11), and Delta(Phis11), which are in turn needed to
    % compute g(x).
    % We will make use of the following chain-rule identities to 
    % compute the derivatives of squared quantities:
    %   1. d(X^2)/dx = 2*X*dX/dx
    %   2. d^2(X^2)/dx^2 = 2[(dX/dx)^2 + X*d^2(X)/dx^2]

    % Fetch phasor quantities and first/second derivatives thereof.
    thetass = lindata.tfThetass(X);
    dthetass = lindata.tfDThetass(1,X);
    d2thetass = lindata.tfDThetass(2,X);
    thetae = lindata.tfThetae(X);
    dthetae = lindata.tfDThetae(1,X);
    d2thetae = lindata.tfDThetae(2,X);
    eta = lindata.tfEta(X);
    deta = lindata.tfDEta(1,X);
    d2eta = lindata.tfDEta(2,X);

    % Compute squared quantities and first/second derivatives thereof.
    thetassSQ = thetass.^2;
    thetaeSQ = thetae.^2;
    etaSQ = eta.^2;
    dthetassSQ = 2*thetass.*dthetass;
    dthetaeSQ = 2*thetae.*dthetae;
    detaSQ = 2*eta.*deta;
    d2thetassSQ = 2*thetass.*d2thetass + 2*dthetass.^2;
    d2thetaeSQ = 2*thetae.*d2thetae + 2*dthetae.^2;
    d2etaSQ = 2*eta.*d2eta + 2*deta.^2;

    % Compute Delta(If), Delta(Idl), Delta(Phis), Delta(Phie),
    % Delta(Eta).
    DeltaIdl = zeros(meta.ns,length(X));
    DeltaIf = zeros(meta.ns,length(X));
    DeltaIfdl = zeros(meta.ns,length(X));
    DeltaPhis = zeros(meta.ns,length(X));
    DeltaPhie = zeros(meta.ns,length(X));
    DeltaEta = zeros(meta.ns,length(X));
    if any(bins.n)
        DeltaIdl(:,bins.n) = -kD*T*thetaeSQ(:,bins.n)./Zdl2n/4;
        DeltaIf(:,bins.n) = 0 ...
            - kD*T*thetaeSQ(:,bins.n)/4/Rct1n ...
            + F*etaSQ(:,bins.n)*Rct2invn/4/R/T;
        DeltaIfdl(:,bins.n) = DeltaIf(:,bins.n) + DeltaIdl(:,bins.n);
        DeltaPhis(:,bins.n) = Rfn*DeltaIfdl(:,bins.n);
        DeltaPhie(:,bins.n) = kD*T*thetaeSQ(:,bins.n)/4;
        DeltaEta(:,bins.n) = -kD*T*thetaeSQ(:,bins.n)/4;
    end
    if any(bins.d|bins.s)
        DeltaIdl(:,bins.d|bins.s) = NaN;
        DeltaIf(:,bins.d|bins.s) = NaN;
        DeltaIfdl(:,bins.d|bins.s) = NaN;
        DeltaPhis(:,bins.d|bins.s) = NaN;
        DeltaPhie(:,bins.d|bins.s) = NaN;
        DeltaEta(:,bins.d|bins.s) = NaN;
    end
    if any(bins.p)
        DeltaIdl(:,bins.p) = -kD*T*thetaeSQ(:,bins.p)./Zdl2p/4;
        DeltaIf(:,bins.p) = -d2Uocpp*thetassSQ(:,bins.p)./Rct1p/4 ...
            - kD*T*thetaeSQ(:,bins.p)/4/Rct1p ...
            + F*etaSQ(:,bins.p)*Rct2invp/4/R/T;
        DeltaIfdl(:,bins.p) = DeltaIf(:,bins.p) + DeltaIdl(:,bins.p);
        DeltaPhis(:,bins.p) = Rfp*DeltaIfdl(:,bins.p);
        DeltaPhie(:,bins.p) = kD*T*thetaeSQ(:,bins.p)/4;
        DeltaEta(:,bins.p) = -d2Uocpp*thetassSQ(:,bins.p)./4 ... 
                           - kD*T*thetaeSQ(:,bins.p)./4;
    end

    % Compute first partials of Delta(If), Delta(Idl), Delta(Phis)
    % w/r/t x(tilde).
    dDeltaIdl = zeros(meta.ns,length(X));
    dDeltaIf = zeros(meta.ns,length(X));
    dDeltaIfdl = zeros(meta.ns,length(X));
    dDeltaPhis = zeros(meta.ns,length(X));
    if any(bins.n)
        dDeltaIdl(:,bins.n) = -kD*T*dthetaeSQ(:,bins.n)./Zdl2p/4;
        dDeltaIf(:,bins.n) = 0 ...
            - kD*T*dthetaeSQ(:,bins.n)/4/Rct1n ...
            + F*detaSQ(bins.n)*Rct2invn/4/R/T;
        dDeltaIfdl(:,bins.n) = dDeltaIf(:,bins.n) + dDeltaIdl(:,bins.n);
        dDeltaPhis(:,bins.n) = Rfn*dDeltaIfdl(:,bins.n);
    end
    if any(bins.d|bins.s)
        dDeltaIdl(:,bins.d|bins.s) = NaN;
        dDeltaIf(:,bins.d|bins.s) = NaN;
        dDeltaIfdl(:,bins.d|bins.s) = NaN;
        dDeltaPhis(:,bins.d|bins.s) = NaN;
    end
    if any(bins.p)
        dDeltaIdl(:,bins.p) = -kD*T*dthetaeSQ(:,bins.p)./Zdl2p/4;
        dDeltaIf(:,bins.p) = -d2Uocpp*dthetassSQ(:,bins.p)./Rct1p/4 ...
                    - kD*T*dthetaeSQ(:,bins.p)/4/Rct1p ...
                    + F*detaSQ(:,bins.p)*Rct2invp/4/R/T;
        dDeltaIfdl(:,bins.p) = dDeltaIf(:,bins.p) + dDeltaIdl(:,bins.p);
        dDeltaPhis(:,bins.p) = Rfp*dDeltaIfdl(:,bins.p);
    end

    % Compute second partials of Delta(If), Delta(Idl), Delta(Phis)
    % w/r/t x(tilde).
    d2DeltaIdl = zeros(meta.ns,length(X));
    d2DeltaIf = zeros(meta.ns,length(X));
    d2DeltaIfdl = zeros(meta.ns,length(X));
    d2DeltaPhis = zeros(meta.ns,length(X));
    if any(bins.n)
        d2DeltaIdl(:,bins.n) = -kD*T*d2thetaeSQ(:,bins.n)./Zdl2p/4;
        d2DeltaIf(:,bins.n) = 0 ...
            - kD*T*d2thetaeSQ(:,bins.n)/4/Rct1n ...
            + F*d2etaSQ(:,bins.n)*Rct2invn/4/R/T;
        d2DeltaIfdl(:,bins.n) = d2DeltaIf(:,bins.n) + d2DeltaIdl(:,bins.n);
        d2DeltaPhis(:,bins.n) = Rfn*d2DeltaIfdl(:,bins.n);
    end
    if any(bins.d|bins.s)
        dDeltaIdl(:,bins.d|bins.s) = NaN;
        dDeltaIf(:,bins.d|bins.s) = NaN;
        dDeltaIfdl(:,bins.d|bins.s) = NaN;
        dDeltaPhis(:,bins.d|bins.s) = NaN;
    end
    if any(bins.p)
        d2DeltaIdl(:,bins.p) = -kD*T*d2thetaeSQ(:,bins.p)./Zdl2p/4;
        d2DeltaIf(:,bins.p) = -d2Uocpp*d2thetassSQ(:,bins.p)./Rct1p/4 ...
                    - kD*T*d2thetaeSQ(:,bins.p)/4/Rct1p ...
                    + F*d2etaSQ(:,bins.p)*Rct2invp/4/R/T;
        d2DeltaIfdl(:,bins.p) = d2DeltaIf(:,bins.p) + d2DeltaIdl(:,bins.p);
        d2DeltaPhis(:,bins.p) = Rfp*d2DeltaIfdl(:,bins.p);
    end

    % Compute g(x(tilde)) and the first and second partials thereof w/r/t
    % x(tilde).
    g = zeros(meta.ns,length(X));
    dg = zeros(meta.ns,length(X));
    d2g = zeros(meta.ns,length(X));
    if any(bins.n|bins.d|bins.s)
        g(:,bins.n|bins.d|bins.s) = NaN;
        dg(:,bins.n|bins.d|bins.s) = NaN;
        d2g(:,bins.n|bins.d|bins.s) = NaN;
    end
    if any(bins.p)
        Z1p = (Zdl2p.*Zs2p)./(Zs2p+Rct1p+Zdl2p)./Zse2p;
        g(:,bins.p) = Z1p.*DeltaIf(:,bins.p)./kappap./psi./T ...
            - DeltaIfdl(:,bins.p)./kappap./psi./T ...
            + DeltaPhis(:,bins.p)./kappap./psi./T./Zse2p;
        dg(:,bins.p) = Z1p.*dDeltaIf(:,bins.p)./kappap./psi./T ...
            - dDeltaIfdl(:,bins.p)./kappap./psi./T ...
            + dDeltaPhis(:,bins.p)./kappap./psi./T./Zse2p;
        d2g(:,bins.p) = Z1p.*d2DeltaIf(:,bins.p)./kappap./psi./T ...
            - d2DeltaIfdl(:,bins.p)./kappap./psi./T ...
            + d2DeltaPhis(:,bins.p)./kappap./psi./T./Zse2p;
    end
    
    data.DeltaIdl = DeltaIdl;
    data.DeltaIf = DeltaIf;
    data.DeltaIfdl = DeltaIfdl;
    data.DeltaPhis = DeltaPhis;
    data.DeltaPhie = DeltaPhie;
    data.DeltaPhise = DeltaPhis - DeltaPhie;
    data.DeltaEta = DeltaEta;
    data.dDeltaIdl = dDeltaIdl;
    data.dDeltaIf = dDeltaIf;
    data.dDeltaIfdl = dDeltaIfdl;
    data.dDeltaPhis = dDeltaPhis;
    data.d2DeltaIdl = d2DeltaIdl;
    data.d2DeltaIf = d2DeltaIf;
    data.d2DeltaIfdl = d2DeltaIfdl;
    data.d2DeltaPhis = d2DeltaPhis;
    data.g = g;
    data.dg = dg;
    data.d2g = d2g;
end