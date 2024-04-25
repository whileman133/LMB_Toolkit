function data = tfLMB(S, model, varargin)
%TFLMB Compute the transfer-functions for a full LMB cell. Linear and
%  second-harmonic response. Faster than using tfXX functions for computing
%  cell impedance, but doesn't implement low/high frequency gains.
%
% -- Usage --
% data = TFLMB(S,model) constructs functions for evalulating the LMB
%   transfer functions at the frequency points specified by the vector S.
%
% data = TFLMB(...,'TdegC',TdegC) performs the computation at temperature 
%    TdegC instead of the default 25degC.
%
% data = TFLMB(...,'Calc11',false) disables calculation of the linear
%   response; only the parameters evalulated at the setpoint are returned.
%
% data = TFLMB(...,'Calc22',true) also calculates the second-harmonic
%   response. The linear response will also be computed, as it is needed to
%   find the second-harmonic response.
%
% -- Changelog --
% 2023.08.31 | Use Warburg params | Wes H.
% 2023.07.02 | Updated for gen2 | Wesley Hileman <whileman@uccs.edu>

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

if isCellModel(model)
    % Ensure we're using legacy LPM.
    model = convertCellModel(model,'LLPM');
end

if isempty(param)
    param = getParameterValues(model,TdegC,socPct,S,calc11,calc22);
end
if calc11 || calc22
    data.h11 = tfLMB11(param,S,param.layerReduction);
end
if calc22
    data.h22 = tfLMB22(data.h11,param.layerReduction);
end
data.param = param;

end % tfLMB()

function p = getParameterValues(model,TdegC,socPct,S,calc11,calc22)
%GETPARAMETERVALUES

T = TdegC+273.15;
R = 8.3144598;      % Molar gas constant [J/mol K]
F = 96485.3329;     % Faraday constant [C/mol]
f = F/R/T;

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
p.nFp         = getParam('pos','nF');
p.tauFp       = getParam('pos','tauF');
p.Rfp         = getParam('pos','Rf');
p.Rdlp        = getParam('pos','Rdl');
p.Cdlp        = getParam('pos','Cdl');
p.nDLp        = getParam('pos','nDL');
p.tauDLp      = getParam('pos','tauDL');
p.theta0p     = getParam('pos','theta0');
p.theta100p   = getParam('pos','theta100');
if isParam('pos','k0')
    p.k0p         = getParam('pos','k0');
end
if isParam('pos','alpha')
    p.alphap      = getParam('pos','alpha');
end
p.Rfn         = getParam('neg','Rf');
p.Rdln        = getParam('neg','Rdl');
p.Cdln        = getParam('neg','Cdl');
p.nDLn        = getParam('neg','nDL');
p.tauDLn      = getParam('neg','tauDL');
p.k0n         = getParam('neg','k0');
p.alphan      = getParam('neg','alpha');
% Code later in this file implements Warburg CPE factor for inert layers;
% not presently used, so we set all nE=1.
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
% Legcacy code below expects qe, kD electrolyte parameters; calculate them
% from tau, W.
p.qep = p.psi*p.T*p.kappap*p.taup/3600;
p.qes = p.psi*p.T*p.kappas*p.taus/3600;
p.qed = p.psi*p.T*p.kappad*p.taud/3600;
p.kD = -p.psi*p.W;
if isnan(p.qed)
    p.qed = 0; % dll combined with sep
end

% Get Uocp, dUocp, and Rct for positive electrode.
if isParam('pos','X') && isParam('pos','U0') && isParam('pos','omega')
    % MSMR model.
    if ~calc22 && ... 
       isParam('pos','Uocp')&&isParam('pos','dUocp')&&isParam('pos','Rct')
        % Linear quantities Uocp, dUocp, and Rct have been pre-computed for speed.
        p.Uocpp = getParam('pos','Uocp');
        p.dUocpp = getParam('pos','dUocp');
        p.Rctp = getParam('pos','Rct');
    elseif isParam('pos','Uocp') && isParam('pos','dUocp') && ... 
           isParam('pos','d2Uocp') && isParam('pos','Rct') && ... 
           isParam('pos','Rct2Inv')
        % Second-harmonic quantities have been precomputed for speed.
        p.Uocpp = getParam('pos','Uocp');
        p.dUocpp = getParam('pos','dUocp');
        p.d2Uocpp = getParam('pos','d2Uocp');
        p.Rctp = getParam('pos','Rct');
        p.Rct2invp = getParam('pos','Rct2Inv');
    else
        % Not all variables have been precomputed; compute here.
        electrode = MSMR(getReg('pos'));
        dataCt = electrode.Rct(getReg('pos'),'npoints',10000,'TdegC',TdegC);
        RctVect = dataCt.Rct;
        UocpVect = dataCt.Uocp;
        dUocpVect = dataCt.dUocp;
        thetaVect = dataCt.theta;
        Rct2InvVect = dataCt.Rct2Inv;
        if ~isParam('pos','Uocp')
            p.Uocpp = interp1(thetaVect,UocpVect,p.thetap,'linear','extrap');
        else
            p.Uocpp = getParam('pos','Uocp');
        end
        if ~isParam('pos','dUocp')
            p.dUocpp = interp1(thetaVect,dUocpVect,p.thetap,'linear','extrap');
        else
            p.dUocpp = getParam('pos','dUocp');
        end
        if ~isParam('pos','Rct')
            p.Rctp = interp1(thetaVect,RctVect,p.thetap,'linear','extrap');
        else
            p.Rctp = getParam('pos','Rct');
        end
        if ~isParam('pos','d2Uocp')
            p.d2Uocpp = interp1(thetaVect,dataCt.d2Uocp,p.thetap,'linear','extrap');
        else
            p.d2Uocpp = getParam('pos','d2Uocp');
        end
        if ~isParam('pos','Rct2inv')
            p.Rct2invp = interp1(thetaVect,Rct2InvVect,p.thetap,'linear','extrap');
        else
            p.Rct2invp = getParam('pos','Rct2inv');
        end
    end
elseif isParam('pos','Uocp') && isParam('pos','dUocp')
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
    % Ds already pre-computed for speed.
    p.Dsp = getParam('pos','Ds');
elseif isParam('pos','Dsref')
    % Ds not precomputed; compute here using Baker-Verbrugge Relationship.
    Dsref = getParam('pos','Dsref');
    p.Dsp = Dsref.*(-f.*thetap.*(1-thetap).*p.dUocpp);
else
    % Ds not precomputed; compute here using lookup table.
    electrode = MSMR(getReg('pos'));
    dsData = electrode.Ds(getReg('pos'),'npoints',10000,'TdegC',TdegC);
    p.Dsp = interp1(dsData.theta,dsData.Ds,p.thetap,'linear','extrap');
end

% Quantities needed to evalulate dc gain.
p.Dseffp = p.Dsp^p.nFp*p.tauFp^(p.nFp-1);
p.Cdleffp = p.Cdlp^p.nDLp*p.tauDLp^(1-p.nDLp); % double-layer cap @ dc
p.Cdleffn = p.Cdln^p.nDLn*p.tauDLn^(1-p.nDLn); % double-layer cap @ dc
p.Csp = -3600*p.Q/p.dUocpp/abs(p.theta100p-p.theta0p); % cell cap @ dc
p.R0p = p.Rctp + 1/15/p.Csp/p.Dseffp;
p.R2p = 1/6/p.sigmap - ... 
        (1+p.W)/3/p.kappap - ... 
        ((p.Rdlp+p.Rfp)*p.Cdleffp^2+(p.R0p+p.Rfp)*p.Csp^2+...
          2*p.Rfp*p.Cdleffp*p.Csp+...
          p.tauDLp*(1-p.nDLp)*p.Cdleffp+...
          p.tauFp*(1-p.nFp)*p.Csp)...
        /(p.Cdleffp+p.Csp)^2;

% Quantities needed to evalulate HF gain.
p.Rseinfn = p.Rfn + p.Rdln*p.Rctn/(p.Rdln+p.Rctn);
p.Rseinfp = p.Rfp + p.Rdlp*p.Rctp/(p.Rdlp+p.Rctp);
p.zetap = sqrt((1/p.sigmap+1/p.kappap)/p.Rseinfp);

% Compute Zse.
ind0 = S==0;  % logical indicies to zero frequencies
p.Zs0 = abs(p.theta100p-p.theta0p)*p.dUocpp/10800/p.Dseffp/p.Q;
p.Dskernp = (p.Dsp./S).*((1+S*p.tauFp)/p.tauFp/p.Dsp).^(1-p.nFp);
p.betap = sqrt(1./p.Dskernp);
p.Zsp = ...
   p.Zs0*( ...
     (p.betap.^2 + 3*(1-p.betap.*coth(p.betap))) ...
       ./(p.betap.^2.*(1-p.betap.*coth(p.betap))) ...
     - 3*p.Dskernp ...
   );
p.Zsp(ind0) = Inf;  % replace NaN at zero frequency
p.Zdlp = p.Rdlp + ((p.Cdlp/p.tauDLp)*(1+p.tauDLp*S)).^(1-p.nDLp)./p.Cdlp./S;
p.Zsep = p.Rfp + 1./(1./(p.Rctp + p.Zsp) + 1./p.Zdlp);
p.Zdln = p.Rdln + ((p.Cdln/p.tauDLn)*(1+p.tauDLn*S)).^(1-p.nDLn)./p.Cdln./S;
p.Zsen = p.Rfn + p.Zdln.*p.Rctn./(p.Zdln+p.Rctn);

% Compute Zse2 (double the frequency!)
S2 = 2*S;
p.Dskern2p = (p.Dsp./S2).*((1+S2*p.tauFp)/p.tauFp/p.Dsp).^(1-p.nFp);
p.beta2p = sqrt(1./p.Dskern2p);
p.Zs2p = ...
   p.Zs0*( ...
     (p.beta2p.^2 + 3*(1-p.beta2p.*coth(p.beta2p))) ...
       ./(p.beta2p.^2.*(1-p.beta2p.*coth(p.beta2p))) ...
     - 3*p.Dskern2p ...
   );
p.Zs2p(ind0) = Inf;  % replace NaN at zero frequency
p.Zdl2p = p.Rdlp + ((p.Cdlp/p.tauDLp)*(1+p.tauDLp*S2)).^(1-p.nDLp)./p.Cdlp./S2;
p.Zse2p = p.Rfp + 1./(1./(p.Rctp + p.Zs2p) + 1./p.Zdl2p);
p.Zdl2n = p.Rdln + ((p.Cdln/p.tauDLn)*(1+p.tauDLn*S2)).^(1-p.nDLn)./p.Cdln./S2;
p.Zse2n = p.Rfn + p.Zdl2n.*p.Rctn./(p.Zdl2n+p.Rctn);

% Additional.
p.split = splitFactory();
p.xlim = xlimFactory();
p.layerReduction = layerReduction;

    function fcn = splitFactory()
        function [Z, bins, regNames] = split1(X)
            X = X(:).';
            bins.n = X==0;
            bins.d = 0<X&X<=1;
            bins.s = 1<X&X<2;
            bins.p = 2<=X&X<=3;  % Some quantities not defined in s, so ensure x=2 is in p!
            Z.n = X(bins.n);
            Z.d = X(bins.d);
            Z.s = X(bins.s)-1;
            Z.p = 3-X(bins.p);
            regNames = cell(1,length(X));
            regNames(bins.n) = {'neg'};
            regNames(bins.d) = {'dll'};
            regNames(bins.s) = {'sep'};
            regNames(bins.p) = {'pos'};
        end
        function [Z, bins, regNames] = split2(X)
            X = X(:).';
            bins.n = X==0;
            bins.d = false(size(X));
            bins.s = 0<X&X<1;
            bins.p = 1<=X&X<=2;  % Some quantities not defined in s, so ensure x=1 is in p!
            Z.n = X(bins.n);
            Z.d = [];
            Z.s = X(bins.s);
            Z.p = 2-X(bins.p);
            regNames = cell(1,length(X));
            regNames(bins.n) = {'neg'};
            regNames(bins.s) = {'sep'};
            regNames(bins.p) = {'pos'};
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
    data.tfIf = @(X)getIf(meta,X,0);
    data.tfPhise = @(X)getPhise(meta,X,0);
    data.tfDPhise = @(r,X)getPhise(meta,X,r);
    data.tfPhiseStar = @(X)getPhiseStar(meta,X,0);
    data.tfThetass = @(X)getThetass(meta,X,0);
    data.tfDThetass = @(r,X)getThetass(meta,X,r);
    data.tfThetassStar = @(X)getThetassStar(meta,X,0);
    data.tfEta = @(X)getEta(meta,X,0);
    data.tfDEta = @(r,X)getEta(meta,X,r);
    data.tfPhie = @(X)getPhie(meta,X);
    data.tfPhieTilde = @(X)getPhieTilde(meta,X);
    data.tfPhisTilde = @(X)getPhisTilde(meta,X);
    data.tfZcell = @()getZcell(meta);
    data.tfZcellStar = @()getZcellStar(meta);
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
            if S(i)==0
                % dc frequency; do not evalulate by solving system of
                % equations; instead, we will compute the dc gain
                % separately in the tfXX functions.
                continue;
            end

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
            if S(i)==0
                % dc frequency; do not evalulate by solving system of
                % equations; instead, we will compute the dc gain
                % separately in the tfXX functions.
                continue;
            end

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

function [Thetae, data] = getThetae(meta,X,r)
    [Z, bins, regNames] = meta.param.split(X);
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

    % Compute dc gain.
    dcGain = zeros(1,length(X));
    p = meta.param;
    psiT = p.psi*p.T;
    Thetae0 = p.taud/2+(1/2+p.kappas/p.kappad)*p.taus+(1/3+p.kappap/p.kappad+p.kappap/p.kappas)*p.taup;
    Thetae0 = Thetae0/(p.kappad*p.taud+p.kappas*p.taus+p.kappap*p.taup)/psiT;
    if any(bins.n)
        dcGain(bins.n) = Thetae0;
    end
    if any(bins.d)
        dcGain(bins.d) = Thetae0-Z.d/p.kappad/psiT;
    end
    if any(bins.s)
        dcGain(bins.s) = Thetae0-(1/p.kappad+Z.s/p.kappas)/psiT;
    end
    if any(bins.p)
        dcGain(bins.p) = Thetae0-(1/p.kappad+1/p.kappas+(1-Z.p).*(1+Z.p)/2/p.kappap)/psiT;
    end
    % Assign dc gain to zero frequencies.
    ind0 = meta.S==0; % logical indicies to dc frequencies
    if any(ind0) && r==0
        Thetae(ind0,:) = dcGain.*ones(sum(ind0),1);
    end % any ind0

    data.dcGain = dcGain;
    data.hfGain = zeros(1,length(X));
    data.res0 = zeros(1,length(X));  % no integrator
    data.regNames = regNames;
end

function [Ifdl, data] = getIfdl(meta,X,r)
    [Z, bins, regNames] = meta.param.split(X);
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

    % Compute dc gain.
    dcGain = zeros(1,length(X));
    if any(bins.n)
        dcGain(bins.n) = 1;
    end
    if any(bins.d|bins.s)
        dcGain(bins.d|bins.s) = NaN;
    end
    if any(bins.p)
        dcGain(bins.p) = -1;
    end
    % Assign dc gain to zero frequencies.
    ind0 = meta.S==0; % logical indicies to dc frequencies
    if any(ind0) && r==0
        Ifdl(ind0,:) = dcGain.*ones(sum(ind0),1);
    end % any ind0

    % Compute hf gain.
    hfGain = zeros(1,length(X));
    if any(bins.n)
        hfGain(bins.n) = 1;
    end
    if any(bins.d|bins.s)
        hfGain(bins.d|bins.s) = NaN;
    end
    if any(bins.p)
        zeta = meta.param.zetap;
        sp = meta.param.sigmap;
        kp = meta.param.kappap;
        Rse = meta.param.Rseinfp;
        hfGain(bins.p) = ...
            -(cosh(zeta*(Z.p-1))/sp + cosh(zeta*Z.p)/kp)/zeta/sinh(zeta)/Rse;
    end

    data.dcGain = dcGain;
    data.hfGain = hfGain;
    data.res0 = zeros(1,length(X));  % no integrator
    data.regNames = regNames;
end

function [If, data] = getIf(meta,X,r)
    [Z, bins, regNames] = meta.param.split(X);
    If = zeros(meta.ns,length(X));

    if any(bins.n)
        if r == 0
            If(:,bins.n) = meta.divFactorIfn;
        else
            If(:,bins.n) = NaN;
        end
    end
    if any(bins.d)
        If(:,bins.d) = NaN;
    end
    if any(bins.s)
        If(:,bins.s) = NaN;
    end
    if any(bins.p)
        If(:,bins.p) = meta.divFactorIfp.*getIfdl(meta,X(bins.p),r);
    end

    % Compute dc gain.
    dcGain = zeros(1,length(X));
    if any(bins.n)
        dcGain(bins.n) = 1;
    end
    if any(bins.d|bins.s)
        dcGain(bins.d|bins.s) = NaN;
    end
    if any(bins.p)
        p = meta.param;
        Cs = p.Csp;
        Cdl = p.Cdleffp;
        dcGain(bins.p) = -Cs/(Cs+Cdl);
    end
    % Assign dc gain to zero frequencies.
    ind0 = meta.S==0; % logical indicies to dc frequencies
    if any(ind0) && r==0
        If(ind0,:) = dcGain.*ones(sum(ind0),1);
    end % any ind0

    % Compute hf gain.
    hfGain = zeros(1,length(X));
    if any(bins.n)
        p = meta.param;
        hfGain(bins.n) = p.Rdln/(p.Rdln+p.Rctn);
    end
    if any(bins.d|bins.s)
        hfGain(bins.d|bins.s) = NaN;
    end
    if any(bins.p)
        zeta = meta.param.zetap;
        sp = meta.param.sigmap;
        kp = meta.param.kappap;
        Rse = meta.param.Rseinfp;
        Rdl = meta.param.Rdlp;
        Rct = meta.param.Rctp;
        hfGain(bins.p) = ...
            -(Rdl/(Rdl+Rct))*(cosh(zeta*(Z.p-1))/sp + cosh(zeta*Z.p)/kp)/zeta/sinh(zeta)/Rse;
    end

    data.dcGain = dcGain;
    data.hfGain = hfGain;
    data.res0 = zeros(1,length(X));  % no integrator
    data.regNames = regNames;
end

function [Phise, data] = getPhiseStar(meta,X,r)
    [Z, bins, regNames] = meta.param.split(X);
    Phise = zeros(meta.ns,length(X));

    % Compute integrator residue.
    p = meta.param;
    Cs = p.Csp;
    Cdl = p.Cdleffp;
    res0p = -1/(Cs+Cdl);
    res0 = zeros(1,length(X));
    if any(bins.p)
        res0(bins.p) = res0p;
    end

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
        Phise(:,bins.p) = meta.param.Zsep.*getIfdl(meta,X(bins.p),r) - res0p./meta.S;
    end

    % Compute dc gain.
    dcGain = zeros(1,length(X));
    if any(bins.n)
        dcGain(bins.n) = meta.param.Rctn+meta.param.Rfn;
    end
    if any(bins.d|bins.s)
        dcGain(bins.d|bins.s) = NaN;
    end
    if any(bins.p)
        sp = p.sigmap;
        kp = p.kappap;
        W = p.W;
        dcGain(bins.p) = p.R2p + ... 
            (2/sp+3/kp-(1/sp+1/kp).*(5-Z.p)./2-W.*(-1-Z.p)/2/kp).*(1-Z.p);
    end
    % Assign dc gain to zero frequencies.
    ind0 = meta.S==0; % logical indicies to dc frequencies
    if any(ind0) && r==0
        Phise(ind0,:) = dcGain.*ones(sum(ind0),1);
    end % any ind0

    % Compute high-frequency gain.
    hfGain = zeros(1,length(X));
    if any(bins.n)
        hfGain(bins.n) = meta.param.Rseinfn;
    end
    if any(bins.d|bins.s)
        hfGain(bins.d|bins.s) = NaN;
    end
    if any(bins.p)
        zeta = meta.param.zetap;
        sp = meta.param.sigmap;
        kp = meta.param.kappap;
        hfGain(bins.p) = ...
            -(cosh(zeta*(Z.p-1))/sp + cosh(zeta*Z.p)/kp)/zeta/sinh(zeta);
    end

    data.res0 = res0;
    data.dcGain = dcGain;
    data.hfGain = hfGain;
    data.regNames = regNames;
end

function [Phise, data] = getPhise(meta,X,r)
    [~, bins, regNames] = meta.param.split(X);
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

    data.regNames = regNames;
end

function [Thetass, data] = getThetassStar(meta,X,r)
%GETTHETASSTILDE Thetass TF without integrator.
    [Z, bins, regNames] = meta.param.split(X);
    Thetass = zeros(meta.ns,length(X));

    % Compute integrator residue.
    p = meta.param;
    Cs = p.Csp;
    Cdldc = p.Cdleffp;
    dUocp = p.dUocpp;
    res0p = -1/dUocp/(Cs+Cdldc);
    res0 = zeros(1,length(X));
    if any(bins.p)
        res0(bins.p) = res0p;
    end

    if any(bins.n|bins.d|bins.s)
        Thetass(:,bins.n|bins.d|bins.s) = NaN;
    end
    if any(bins.p)
        If = getIf(meta,X(bins.p),r);
        Thetass(:,bins.p) = (p.Zsp./dUocp).*If - res0p./meta.S;
    end

    % Compute dc gain.
    dcGain = zeros(1,length(X));
    p = meta.param;
    if any(bins.n|bins.d|bins.s)
        dcGain(bins.n|bins.d|bins.s) = NaN;
    end
    if any(bins.p)
        sp = p.sigmap;
        kp = p.kappap;
        W = p.W;
        dUocp = p.dUocpp;
        Cdldc = p.Cdleffp;
        Cs = p.Csp;
        Rct = p.Rctp;
        Rf = p.Rfp;
        R2 = p.R2p;
        R3 = (2/sp+3/kp-(1/sp+1/kp).*(5-Z.p)./2-W.*(-1-Z.p)/2/kp).*(1-Z.p);
        dcGain(bins.p) = ( ... 
           R2+R3+Rf + (Cs/(Cs+Cdldc))*Rct ...
        )/dUocp;
    end
    % Assign dc gain to zero frequencies.
    ind0 = meta.S==0; % logical indicies to dc frequencies
    if any(ind0) && r==0
        Thetass(ind0,:) = dcGain.*ones(sum(ind0),1);
    end % any ind0

    % HF gain.
    hfGain = zeros(1,length(X));

    data.dcGain = dcGain;
    data.hfGain = hfGain;
    data.res0 = res0;
    data.regNames = regNames;
end

function [Thetass, data] = getThetass(meta,X,r)
%GETTHETASS Thetass TF with integrator.
    [~, bins, regNames] = meta.param.split(X);
    Thetass = zeros(meta.ns,length(X));

    if any(bins.n|bins.d|bins.s)
        Thetass(:,bins.n|bins.d|bins.s) = NaN;
    end

    if any(bins.p)
        p = meta.param;
        If = meta.divFactorIfp.*getIfdl(meta,X(bins.p),r);
        Thetass(:,bins.p) = (p.Zsp./p.dUocpp).*(If  );
    end

    data.regNames = regNames;
end

function [Eta, data] = getEta(meta,X,r)
    [Z, bins, regNames] = meta.param.split(X);
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

     % Compute dc gain.
    dcGain = zeros(1,length(X));
    if any(bins.n)
        dcGain(bins.n) = meta.param.Rctn;
    end
    if any(bins.d|bins.s)
        dcGain(bins.d|bins.s) = NaN;
    end
    if any(bins.p)
        p = meta.param;
        Cs = p.Csp;
        Cdl = p.Cdleffp;
        Rct = p.Rctp;
        dcGain(bins.p) = -Rct*Cs/(Cs+Cdl);
    end
    % Assign dc gain to zero frequencies.
    ind0 = meta.S==0; % logical indicies to dc frequencies
    if any(ind0) && r==0
        Eta(ind0,:) = dcGain.*ones(sum(ind0),1);
    end % any ind0

    % Compute high-frequency gain.
    hfGain = zeros(1,length(X));
    if any(bins.n)
        hfGain(bins.n) = meta.param.Rseinfn-meta.param.Rfn;
    end
    if any(bins.d|bins.s)
        hfGain(bins.d|bins.s) = NaN;
    end
    if any(bins.p)
        zeta = meta.param.zetap;
        sp = meta.param.sigmap;
        kp = meta.param.kappap;
        Rf = meta.param.Rfp;
        Rseinf = meta.param.Rseinfp;
        hfGain(bins.p) = ...
            -(1-Rf/Rseinf)*( ...
                cosh(zeta*(Z.p-1))/sp + cosh(zeta*Z.p)/kp ...
            )/zeta/sinh(zeta);
    end

    data.dcGain = dcGain;
    data.hfGain = hfGain;
    data.res0 = zeros(1,length(X));  % no integrator
    data.regNames = regNames;
end

function [Phie, data] = getPhie(meta,X)
    PhieTilde = getPhieTilde(meta,X);
    Phise0 = getPhise(meta,0,0);
    Phie = PhieTilde - Phise0;
    data = struct;
end

function [Phie, data] = getPhieTilde(meta,X)
    [Z, bins, regNames] = meta.param.split(X);
    Phie1 = zeros(meta.ns,length(X));

    kd = meta.param.kappad;
    ks = meta.param.kappas;
    kp = meta.param.kappap;
    sp = meta.param.sigmap;
    W = meta.param.W;
    zeta = meta.param.zetap;

    if any(bins.n)
        Phie1(:,bins.n) = 0;
    end
    if any(bins.d)
        Phie1(:,bins.d) = -( Z.d/kd ).*ones(meta.ns,1);
    end
    if any(bins.s)
        Phie1(:,bins.s) = -( 1/kd + Z.s/ks ).*ones(meta.ns,1);
    end
    if any(bins.p)
        L1P = meta.Lambda1P;
        L1Psq = L1P.^2;
        L2P = meta.Lambda2P;
        L2Psq = L2P.^2;
        Phie1(:,bins.p) = -1/kd - 1/ks - (1-Z.p)/kp ...
            - (meta.j1p./L1Psq./kp).*(exp(L1P.*(Z.p-1))-L1P.*(Z.p-1)-1) ...
            - (meta.j2p./L1Psq./kp).*(exp(-L1P.*Z.p)+L1P.*(Z.p-1).*exp(-L1P)-exp(-L1P)) ...
            - (meta.j3p./L2Psq./kp).*(exp(L2P.*(Z.p-1))-L2P.*(Z.p-1)-1) ...
            - (meta.j4p./L2Psq./kp).*(exp(-L2P.*Z.p)+L2P.*(Z.p-1).*exp(-L2P)-exp(-L2P));
    end
    Thetae = getThetae(meta,X,0);
    Thetae0 = getThetae(meta,0,0);
    Phie2 = -meta.param.kD*meta.T*(Thetae-Thetae0);
    Phie = Phie1 + Phie2;
    data.Phie1 = Phie1;
    data.Phie2 = Phie2;

    % Compute dc gain.
    dcGain = zeros(1,length(X));
    if any(bins.n)
        dcGain(bins.n) = 0;
    end
    if any(bins.d)
        dcGain(bins.d) = -(1+W)*Z.d/kd;
    end
    if any(bins.s)
        dcGain(bins.s) = -(1+W)*(1/kd+Z.s/ks);
    end
    if any(bins.p)
        dcGain(bins.p) = -(1+W)*(1/kd+1/ks+(1-Z.p).*(1+Z.p)/2/kp);
    end
    % Assign dc gain to zero frequencies.
    ind0 = meta.S==0; % logical indicies to dc frequencies
    if any(ind0)
        Phie(ind0,:) = dcGain.*ones(sum(ind0),1);
    end % any ind0

    % Compute hf gain.
    hfGain = zeros(1,length(X));
    if any(bins.n)
        hfGain(bins.n) = 0;
    end
    if any(bins.d)
        hfGain(bins.d) = -Z.d/kd;
    end
    if any(bins.s)
        hfGain(bins.s) = -(1/kd + Z.s/ks);
    end
    if any(bins.p)
        hfGain(bins.p) = -( ...
            1/kd + 1/ks + ...
            (1-Z.p)/(kp+sp) + ...
            (sp/(sp+kp))*( ...
              1/sp + cosh(zeta)/kp - ...
              (cosh(zeta*(Z.p-1))/sp + cosh(zeta*Z.p)/kp) ...
            )/zeta/sinh(zeta) ...
       );
    end

    data.dcGain = dcGain;
    data.hfGain = hfGain;
    data.res0 = zeros(1,length(X));  % no integrator
    data.regNames = regNames;
end

function [Phis, data] = getPhisTilde(meta,X)
    [Z, bins, regNames] = meta.param.split(X);
    Phis = zeros(meta.ns,length(X));

    if any(bins.n)
        Phis(:,bins.n) = 0;
    end
    if any(bins.d|bins.s)
        Phis(:,bins.d|bins.s) = NaN;
    end
    if any(bins.p)
        sp = meta.param.sigmap;
        L1P = meta.Lambda1P;
        L2P = meta.Lambda2P;
        j1 = meta.j1p;
        j2 = meta.j2p;
        j3 = meta.j3p;
        j4 = meta.j4p;
        Phis(:,bins.p) = ( ...
            Z.p + ...
            j1.*exp(-L1P).*(exp(L1P.*Z.p)-L1P.*Z.p-1)./L1P.^2 + ...
            j2.*(exp(-L1P.*Z.p)+L1P.*Z.p-1)./L1P.^2 + ...
            j3.*exp(-L2P).*(exp(L2P.*Z.p)-L2P.*Z.p-1)./L2P.^2 + ...
            j4.*(exp(-L2P.*Z.p)+L2P.*Z.p-1)./L2P.^2 ...
        )./sp;
    end

    % Compute dc gain.
    dcGain = zeros(1,length(X));
    if any(bins.n)
        dcGain(bins.n) = 0;
    end
    if any(bins.d|bins.s)
        dcGain(bins.d|bins.s) = NaN;
    end
    if any(bins.p)
        sp = meta.param.sigmap;
        dcGain(bins.p) = -Z.p.*(Z.p-2)/2/sp;
    end
    % Assign dc gain to zero frequencies.
    ind0 = meta.S==0; % logical indicies to dc frequencies
    if any(ind0)
        Phis(ind0,:) = dcGain.*ones(sum(ind0),1);
    end % any ind0

    % Compute hf gain.
    hfGain = zeros(1,length(X));
    if any(bins.n)
        hfGain(bins.n) = 0;
    end
    if any(bins.d|bins.s)
        hfGain(bins.d|bins.s) = NaN;
    end
    if any(bins.p)
        kp = meta.param.kappap;
        sp = meta.param.sigmap;
        zeta = meta.param.zetap;
        hfGain(bins.p) = ( ...
            Z.p/(sp+kp) + ...
            (kp/(sp+kp))*( ...
              1/kp + cosh(zeta)/sp - ...
              (cosh(zeta*(Z.p-1))/sp + cosh(zeta*Z.p)/kp) ...
            )/zeta/sinh(zeta) ...
       );
    end

    data.dcGain = dcGain;
    data.hfGain = hfGain;
    data.res0 = zeros(1,length(X));  % no integrator
    data.regNames = regNames;
end

function [Zcell, data] = getZcell(meta)
    Phise = getPhise(meta,[meta.param.xlim.min meta.param.xlim.max],0);
    Phise0 = Phise(:,1);
    PhiseL = Phise(:,2);
    PhieL = getPhieTilde(meta,meta.param.xlim.max);
    Zcell = -(PhiseL + PhieL - Phise0);
    data.Zneg = Phise0;
    data.Zpos = -PhiseL;
    data.Zel = -PhieL;
end

function [Zcell, data] = getZcellStar(meta)
    [PhiseStar, PhiseStarData] = getPhiseStar(meta,[meta.param.xlim.min meta.param.xlim.max],0);
    PhiseStar0 = PhiseStar(:,1);
    PhiseStarL = PhiseStar(:,2);
    [PhieL, PhieLData] = getPhieTilde(meta,meta.param.xlim.max);
    Zcell = -(PhiseStarL + PhieL - PhiseStar0);
    Rcell0 = -(PhiseStarData.dcGain(2)+PhieLData.dcGain-PhiseStarData.dcGain(1));
    RcellInf = -(PhiseStarData.hfGain(2)+PhieLData.hfGain-PhiseStarData.hfGain(1));
    data.Zneg = PhiseStar0;
    data.Zpos = -PhiseStarL;
    data.Zel = -PhieL;
    data.dcGain = Rcell0;
    data.hfGain = RcellInf;
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
    opts = bvpset('AbsTol',1e-9,'RelTol',1e-4);

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
    data.tfPhieTilde = @(X)getPhie22(meta,X,true);
    data.tfThetass = @(X)getThetass22(meta,X,true);
    data.tfEta = @(X)getEta22(meta,X,true);
    data.tfZcell = @()getZcell22(meta);
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

        % Solve 4th-order BVP in Thetae at each frequency.
        % Allocate storage for solution.
        Thetae22pInit = zeros(ns,length(zsoln));
        Thetae22sInit = zeros(ns,length(zsoln));
        Thetae22p = zeros(ns,length(zsoln));
        Thetae22s = zeros(ns,length(zsoln));
        d1Thetae22p = zeros(ns,length(zsoln));
        d1Thetae22s = zeros(ns,length(zsoln));
        d2Thetae22p = zeros(ns,length(zsoln));  % only available in pos
        for k = 1:ns
            tau1p = tau12p(k);
            tau2p = tau22p(k);
            mu2p = mu22p(k);
    
            % Initialize BVP solver.
            sol = bvpinit(zmesh,@(z)ThetaeInitial([],z));
            Thetae22pInit(k,:) = interp1(sol.x,sol.y(ind.Thetaep,:),zsoln,'linear','extrap');
            Thetae22sInit(k,:) = interp1(sol.x,sol.y(ind.Thetaes,:),zsoln,'linear','extrap');

            % Solve BVP.
            sol = bvp5c(@ThetaeODE, @ThetaeBC, sol, opts);
            Thetae = deval(sol,zsoln,[
                ind.Thetaep
                ind.d1Thetaep
                ind.d2Thetaep
                ind.Thetaes
                ind.d1Thetaes
            ].');
            Thetae22p(k,:) = Thetae(1,:);
            d1Thetae22p(k,:) = Thetae(2,:);
            d2Thetae22p(k,:) = Thetae(3,:);
            Thetae22s(k,:) = Thetae(4,:);
            d1Thetae22s(k,:) = Thetae(5,:);
        end

        soln.xsoln = [zsoln(1:end-1) 1+zsoln(1:end)]; % x=1+z in pos
        soln.Thetae22Init = [Thetae22sInit(:,1:end-1) Thetae22pInit];
        soln.Thetae22 = [Thetae22s(:,1:end-1) Thetae22p];
        soln.d1Thetae22 = [d1Thetae22s(:,1:end-1) d1Thetae22p];
        soln.d2Thetae22 = [nan(ns,length(zsoln)-1) d2Thetae22p];
        soln.i1Thetae22p = cumtrapz(zsoln,Thetae22p.').';
        soln.i2Thetae22p = cumtrapz(zsoln,soln.i1Thetae22p.').';
        soln.i2Thetae22 = [nan(ns,length(zsoln)-1) soln.i2Thetae22p];

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
    
        function yinit = ThetaeInitial(~,Z)
            % Generate initial guess for the Thetae y-vector.

            % Trivial solution. (Does not satisfy some of the pos BCs.)
            yinit = zeros(length(fieldnames(ind)),length(Z));
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
    X = X(:).';
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
            - psi*T*d1ThetaeSP.*(X(bins.p)-meta.param.xlim.sp);
    end

    Phie = Phie1 + Phie2;
    if delinearize
        Phie = Phie + interp1(meta.xsoln,meta.DeltaPhie.',X).';
    end

    data.Phie1 = Phie1;
    data.Phie2 = Phie2;
end

function [Zcell, data] = getZcell22(meta)
    Phise = getPhise22(meta,[meta.param.xlim.min meta.param.xlim.max],true);
    Phise0 = Phise(:,1);
    PhiseL = Phise(:,2);
    PhieL = getPhie22(meta,meta.param.xlim.max,true);
    Zcell = -(PhiseL + PhieL - Phise0);
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