% function [C,Lambda,J,Z,Rct] = tfCommon(s,cellData)
% 
% Inputs:
%   s        = vector of 1j*w where w is vector of frequencies at which
%              transfer functions (TFs) are to be evaluated
%   cellData = data structure containing all parameter values for a cell at
%              a given setpoint, most likely the output of evalSetpoint.m
% Outputs:
%   C        = a matrix of c1n,c2n,c3n,c4n,c1s,c11s,c2s,c1p,c2p,c3p,c4p
%              coefficients required to evaluate TFs 
%   Lambda   = a matrix of the Lambda values required to evaluate TFs
%   J        = a matrix of the "j" values (e.g., required to evaluate ifdl)
%   Z        = a matrix of interface impedances
%   Rct      = a matrix of charge-transfer resistances
%
% This utility function evaluates terms that are common to implementing
% most of the cell transfer functions. It is unlikely that this function
% will ever need to be invoked by the user. It is invoked directly by the
% TF function routines that require its outputs.
%
function [C,Lambda,J,Z,Rct] = tfCommon(s,cellData)
  s = s(:).';  % Force row vector! 

  % first, check to see if we have already computed the relevant data...
  % don't recompute unless necessary!
  if ( ...
      isfield(cellData,'common') && ...
      isfield(cellData.common,'s') && ...
      isequal(cellData.common.s(:),s(:)) ...
  )
    C = cellData.common.C;
    Lambda = cellData.common.L;
    J = cellData.common.J;
    Z = cellData.common.Z;
    Rct = cellData.common.Rct;
    return;
  end

  % Set up negative- and positive-electrode constants
  F = cellData.const.F;
  R = cellData.const.R;
  T = cellData.const.T;
  f = F/(R*T);

  % Set up constant values
  psi = cellData.const.psi;
  kappaD = cellData.const.kD;
  Q = cellData.const.Q;  

  % set up dead-Li layer constants...
  kappaDL = cellData.DL.kappa;
  qeDL = cellData.DL.qe;
  if isfield(cellData.DL,'nE')
    nEdl = cellData.DL.nE;
  else
    nEdl = 1;
  end
  
  % set up separator constants
  kappas = cellData.sep.kappa;
  qes = cellData.sep.qe;
  if isfield(cellData.sep,'nE')
    nEs = cellData.sep.nE;
  else
    nEs = 1;
  end

  % set up positive-electrode constants...
  socp = cellData.pos.soc;
  sigmap = cellData.pos.sigma;
  kappap = cellData.pos.kappa;
  qep = cellData.pos.qe;
  Dsp = cellData.pos.Ds;
  DeltaQp = abs(cellData.pos.theta100 - cellData.pos.theta0);
  csmaxp = 10800*Q*Dsp/DeltaQp;
  nF = cellData.pos.nF;
  Rf = cellData.pos.Rf; 
  Uocpp = cellData.pos.Uocp;
  dUocpp = cellData.pos.dUocp;
  k0p = cellData.pos.k0;
  alphap = cellData.pos.alpha;
  Rdlp = cellData.pos.Rdl;
  Cdlp = cellData.pos.Cdl;
  nDLp = cellData.pos.nDL;  
  wDLp = cellData.pos.wDL;

  % Compute MSMR charge-transfer resistance (independent of frequency!)
  k0 = k0p;
  alpha = alphap; % charge-transfer alpha
  if length(k0) > 1 % Then, this is an MSMR-kinetics model
    Xp = cellData.pos.X;
    U0p = cellData.pos.U0;
    Wp = cellData.pos.omega;

    Uref = Uocpp;
    X    = Xp;
    U0   = U0p;
    W    = Wp;

    % Compute fractions xj and i0 of each gallary (vectors)
    x    = X./(1+exp(f*(Uref-U0)./W)); % sum(xjn) = thetan
    i0   = sum(k0.*((x).^(W.*alpha)).*((X-x).^(W.*(1-alpha))));  
  else % This is a standard Butler-Volmer-kinetics model
    i0   = k0.*((1 - socp)).^(1 - alpha).*socp.^(alpha);   
  end  
  Rctp = R*T/i0/F;

  % Compute positive electrode interface impedance.
  beta = sqrt((s/Dsp).^nF);
  Zsp = (dUocpp/csmaxp)...
       *((beta.^2 + 3*(1-beta.*coth(beta)))...
       ./(beta.^2.*(1-beta.*coth(beta))) - 3*Dsp./s);  % diffusion impedance
  Zdlp = Rdlp + (Cdlp*s+wDLp).^(1-nDLp)./(Cdlp^(2-nDLp)*s);  % double-layer
  Zsep = Rf + 1./(1./(Rctp + Zsp) + 1./Zdlp); % solid-electrolyte impedance
  Isplitp = 1./(1+(Rctp + Zsp)./Zdlp);

  % Pre-allocate storage matricies.
  C = zeros(8,length(s));
  Lambda = zeros(4,length(s));
  j1p = zeros(1,length(s));
  j2p = zeros(1,length(s));
  j3p = zeros(1,length(s));
  j4p = zeros(1,length(s));
  Rctp_ = zeros(1,length(s));

for i=1:length(s)
  % Set up solution variables
  %Separator
  Lambda1s = sqrt(((3600*qes/(psi*T)/kappas)*s(i))^nEs); 
%   if ~isempty(ind0), Lambda1s(ind0)=0; end

  % Dead Li layer
  Lambda1DL = sqrt(((3600*qeDL/(psi*T)/kappaDL)*s(i))^nEdl); 
%   if ~isempty(ind0), Lambda1DL(ind0)=0; end

  mu1p = (1/kappap/(psi*T))./Zsep(i);
  mu2p = (kappaD*T/(psi*T)/kappap)./Zsep(i) - ((3600*qep/(psi*T)/kappap).*s(i));
  tau1p = (1/sigmap + 1/kappap)./Zsep(i) - mu2p; 
  tau2p = ((3600*qep/(psi*T)/kappap)*s(i)) .* (1/sigmap + 1/kappap)./Zsep(i); 
  Lambda1p = sqrt(0.5*(tau1p-sqrt(tau1p.^2-4*tau2p))); 
  Lambda2p = sqrt(0.5*(tau1p+sqrt(tau1p.^2-4*tau2p))); 
  lambda1p = Lambda1p.^3+mu2p.*Lambda1p;
  lambda2p = Lambda2p.^3+mu2p.*Lambda2p;
 
  c1 = [Lambda1DL*exp(-Lambda1DL), -Lambda1DL, 0, 0, 0, 0, 0, 0,];

  c2 = [kappaDL*Lambda1DL, -kappaDL*Lambda1DL*exp(-Lambda1DL),...
    -kappas*Lambda1s*exp(-Lambda1s), kappas*Lambda1s, 0, 0, 0, 0];

  c3 = [1, exp(-Lambda1DL),  -exp(-Lambda1s), -1, 0, 0, 0, 0];

  c4 = [0, 0, kappas*Lambda1s, -kappas*Lambda1s*exp(-Lambda1s), kappap*Lambda1p,...
    -kappap*Lambda1p*exp(-Lambda1p), kappap*Lambda2p, -kappap*Lambda2p*exp(-Lambda2p)];

  c5 = [0, 0, 1, exp(-Lambda1s), -1, -exp(-Lambda1p), -1, -exp(-Lambda2p) ];

  c6 = [ 0, 0, 0, 0, lambda1p, -lambda1p*exp(-Lambda1p), lambda2p, -lambda2p*exp(-Lambda2p)];

  c7 = [ 0, 0, 0, 0, Lambda1p*exp(-Lambda1p), -Lambda1p, Lambda2p*exp(-Lambda2p), -Lambda2p];

  c8 = [ 0, 0, 0, 0,exp(-Lambda1p)*Lambda1p^3, -(Lambda1p)^3, exp(-Lambda2p)*(Lambda2p)^3, -(Lambda2p)^3];

  b = [-1/(kappaDL*(psi*T)); 0; 0; 0; 0; mu1p/kappap; 0; -mu1p/sigmap];

  A = [c1;c2;c3;c4;c5;c6;c7;c8];

%   lastwarn('', '');
  Cs = A\b;
%   [~, warnId] = lastwarn();
%   if ~isempty(warnId)
%       disp(A);
%       disp(s(i)/1j/2/pi);
%   end

  c1p = Cs(5);
  c2p = Cs(6);
  c3p = Cs(7);
  c4p = Cs(8);
  
  j1p(i) = c1p.*(-(psi*T)*kappap*Lambda1p.^2+3600*qep*s(i));
  j2p(i) = c2p.*(-(psi*T)*kappap*Lambda1p.^2+3600*qep*s(i));
  j3p(i) = c3p.*(-(psi*T)*kappap*Lambda2p.^2+3600*qep*s(i));
  j4p(i) = c4p.*(-(psi*T)*kappap*Lambda2p.^2+3600*qep*s(i));
  Lambda(:,i) = [Lambda1DL;Lambda1s;Lambda1p;Lambda2p]; 
  C(:,i) = Cs;
  Rctp_(i) = Rctp;
end

  k0neg = cellData.neg.k0;
  Rctn = R*T/(k0neg*F);

  J = [j1p;j2p;j3p;j4p];
  Z = [Zsep;Zsp;Isplitp]; 
  Rct = [Rctn Rctp_(1)];           
end