% function [C,Lambda,J,Z,Rct,param] = tfCommon(s,cellData)
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
%   param    = flat structure of parameter values required for TF eval.
%
% This utility function evaluates terms that are common to implementing
% most of the cell transfer functions. It is unlikely that this function
% will ever need to be invoked by the user. It is invoked directly by the
% TF function routines that require its outputs.
%
% -- Changelog --
% 2023.09.17 | Use tfLMB() to evalulate parameter values | Wes H.
% 2023.09.14 | Use Warburg parameters | Wesley Hileman <whileman@uccs.edu>
%
function [C,Lambda,J,Z,Rct,param] = tfCommon(s,cellData)
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
        param = cellData.common.param;
        return;
    end
    
    % Fetch cell parameters from tfLMB function. This keeps all of the
    % parameter evalulation in a single place. Modify tfLMB to change how
    % the parameter values are evalulated.
    data = tfLMB(s,cellData, ...
        'Calc11',false,'Calc22',false, ...
        'TdegC',cellData.const.T-273.15,'socPct',cellData.const.soc*100);
    p = data.param;
    T = p.T;
    psi = p.psi;
    kappaD = p.kD;
    kappad = p.kappad;
    qed = p.qed;
    nEdl = p.nEdl;
    kappas = p.kappas;
    qes = p.qes;
    nEs = p.nEs;
    sigmap = p.sigmap;
    kappap = p.kappap;
    qep = p.qep;
    Zsp = p.Zsp(:).';    % ! important, need to be row vectors
    Zdlp = p.Zdlp(:).';
    Zsep = p.Zsep(:).';
    Rctp = p.Rctp;
    Rctn = p.Rctn;
    Isplitp = 1./(1+(Rctp + Zsp)./Zdlp);
    
    % Pre-allocate storage matricies.
    C = zeros(8,length(s));
    Lambda = zeros(4,length(s));
    j1p = zeros(1,length(s));
    j2p = zeros(1,length(s));
    j3p = zeros(1,length(s));
    j4p = zeros(1,length(s));
    
    for i=1:length(s)
        Lambda1s = sqrt(((3600*qes/(psi*T)/kappas)*s(i))^nEs); 
        Lambda1DL = sqrt(((3600*qed/(psi*T)/kappad)*s(i))^nEdl); 
        mu1p = (1/kappap/(psi*T))./Zsep(i);
        mu2p = (kappaD*T/(psi*T)/kappap)./Zsep(i) - ((3600*qep/(psi*T)/kappap).*s(i));
        tau1p = (1/sigmap + 1/kappap)./Zsep(i) - mu2p; 
        tau2p = ((3600*qep/(psi*T)/kappap)*s(i)) .* (1/sigmap + 1/kappap)./Zsep(i); 
        Lambda1p = sqrt(0.5*(tau1p-sqrt(tau1p.^2-4*tau2p))); 
        Lambda2p = sqrt(0.5*(tau1p+sqrt(tau1p.^2-4*tau2p))); 
        lambda1p = Lambda1p.^3+mu2p.*Lambda1p;
        lambda2p = Lambda2p.^3+mu2p.*Lambda2p;
        
        c1 = [Lambda1DL*exp(-Lambda1DL), -Lambda1DL, 0, 0, 0, 0, 0, 0,];
        c2 = [kappad*Lambda1DL, -kappad*Lambda1DL*exp(-Lambda1DL), -kappas*Lambda1s*exp(-Lambda1s), kappas*Lambda1s, 0, 0, 0, 0];
        c3 = [1, exp(-Lambda1DL),  -exp(-Lambda1s), -1, 0, 0, 0, 0];
        c4 = [0, 0, kappas*Lambda1s, -kappas*Lambda1s*exp(-Lambda1s), kappap*Lambda1p, -kappap*Lambda1p*exp(-Lambda1p), kappap*Lambda2p, -kappap*Lambda2p*exp(-Lambda2p)];
        c5 = [0, 0, 1, exp(-Lambda1s), -1, -exp(-Lambda1p), -1, -exp(-Lambda2p) ];
        c6 = [ 0, 0, 0, 0, lambda1p, -lambda1p*exp(-Lambda1p), lambda2p, -lambda2p*exp(-Lambda2p)];
        c7 = [ 0, 0, 0, 0, Lambda1p*exp(-Lambda1p), -Lambda1p, Lambda2p*exp(-Lambda2p), -Lambda2p];
        c8 = [ 0, 0, 0, 0,exp(-Lambda1p)*Lambda1p^3, -(Lambda1p)^3, exp(-Lambda2p)*(Lambda2p)^3, -(Lambda2p)^3];
        b = [-1/(kappad*(psi*T)); 0; 0; 0; 0; mu1p/kappap; 0; -mu1p/sigmap];
        A = [c1;c2;c3;c4;c5;c6;c7;c8];
        Cs = A\b;
        
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
    end % for
    
    J = [j1p;j2p;j3p;j4p];
    Z = [Zsep;Zsp;Isplitp]; 
    Rct = [Rctn Rctp];      
    param = p;
end