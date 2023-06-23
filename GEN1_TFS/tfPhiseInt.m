% function [phiseTF,aux] = tfPhiseInt(s,locs,cellData) 
% 
% Inputs:
%   s        = vector of 1j*w where w is vector of frequencies at which
%              the transfer function (TF) is to be evaluated
%   locs     = vector of locations in normalized coordinates (0..3) at
%              which the TF is to be evaluated
%   cellData = data structure containing all parameter values for a cell at
%              a given setpoint, most likely the output of evalSetpoint.m
% Outputs:
%   phiseTF  = a matrix of the TF evaluated at every combination of "s" and
%              "locs" -- each row holds the TF for one location and all
%              frequencies; there is one row for every location
%   aux      = other outputs sometimes needed by calling routines (e.g.,
%              the xRAs). These outputs include dc-gains, high-frequency
%              gains, integrator residues, a record of the TF name, etc.
%
% This function evaluates the interphase-potential-difference (solid
% potential minus electrolyte potential, WITH integrator) TF. If you wish
% to evaluate the TF without the integrator, use tfPhise.m instead.
%
% Copyright (c) 2021 by Gregory L. Plett and M. Scott Trimboli of the
% University of Colorado Colorado Springs (UCCS). This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 Intl.
% License, v. 1.0. It is provided "as is", without express or implied
% warranty, for educational and informational purposes only.

function [phiseTF,aux] = tfPhiseInt(s,locs,cellData) %#ok<*NASGU>
  % Set up some constants
  F = cellData.const.F;
  R = cellData.const.R;
  T = cellData.const.T;

  % Initialize outputs
  s = s(:).';          % force "s" to be a row vector
  phiseTF = zeros(length(locs),length(s));
  aux.dcGain = Inf(length(locs),1);      % low-frequency gain
  aux.res0 = 0*locs(:);                  % integrator residue
  aux.names = cell(length(locs),1);      % TF names
  aux.xLoc = zeros(1,length(locs));      % TF locations  
  
  % Get common components... (next three lines needed to keep MATLAB happy)
  L1DL=[];L1s=[];L1p=[];L2p=[];
  c1DL=[];c2DL=[];
  c1s=[];c2s=[];
  c1p=[];c2p=[];c3p=[];c4p=[];
  j1p=[];j2p=[];j3p=[];j4p=[];
  Zsep=[];
  Zsp=[];
  Isn=[];Isp=[]; 
  eval(cellData.common.ind); % set c##,L##,J##,Z## indices

  [~,L,J,Z,Rct] = tfCommon(s,cellData); 

  aux.hfGain = calcDterm(locs,Rct(2));      % high-frequency gain of cathode

  indNeg = find(locs == 0);
  if ~isempty(indNeg)
    Rfn = cellData.neg.Rf;
    Rctn = cellData.common.Rct(1);
    Rdln = cellData.neg.Rdl;
    Cdln = cellData.neg.Cdl;
    nDLn = cellData.neg.nDL;  
    wDLn = cellData.neg.wDL;

    F = cellData.const.F;
    R = cellData.const.R;
    T = cellData.const.T;
    f = F/(R*T);
    s = [1e-8 s];

    Zdln = Rdln + (Cdln*s+wDLn).^(1-nDLn)./(Cdln^(2-nDLn)*s);
    Zsfn = (Zdln./(Zdln + Rctn)).*(Rctn);
    s = s(2:end);
    
    phiseTF(indNeg,:) = Zsfn(2:end) + Rfn;
    [aux.names{indNeg}] = deal('negPhise');
    aux.xLoc(indNeg) = locs(indNeg);  
    aux.hfGain(indNeg) = Rctn*Rdln/(Rdln+Rctn) + Rfn;

     inds0 = find(s==0,1);
     if ~isempty(inds0)
        aux.dcGain(indNeg) = Zsfn(1) + Rfn;
        phiseTF(indNeg,inds0) = aux.dcGain(indNeg); 
     end
  end
  
  indPos = find(locs >= 2);
  if ~isempty(indPos)
    for k = 1:length(indPos)
      z = 3-locs(indPos(k));
      ifdl_tf = ...
        J(j1p,:).*exp(L(L1p,:)*(z-1)) + J(j2p,:).*exp(-L(L1p,:)*z) + ...
        J(j3p,:).*exp(L(L2p,:)*(z-1)) + J(j4p,:).*exp(-L(L2p,:)*z);      
      phiseTF(indPos(k),:) = Z(Zsep,:).*ifdl_tf;
      aux.names{indPos(k)} = 'posPhiseInt';
      aux.xLoc(indPos(k)) = 3-z;                                        
    end 
  end

  % Computes high-frequency gains (f->inf)
  function dTerm = calcDterm(locs,Rct)
    dTerm = 0*locs(:); % dummy initialize


    % set up positive-electrode constants...
    Rfp = cellData.pos.Rf;
    Cdlp = cellData.pos.Cdl;     
    Rctp = Rct;
    if Cdlp == 0 % no double layer
      ZseInfp = Rfp + Rctp;
    else
      ZseInfp = Rfp + 1./(1./Rctp + 1./cellData.pos.Rdl);
    end
    sigmap = cellData.pos.sigma;
    kappap = cellData.pos.kappa;

    indPos = find(locs >= 2);
    nuInfp = sqrt((1/sigmap+1/kappap)/ZseInfp);
    if ~isempty(indPos)
      for kk = 1:length(indPos)
        z = 3-locs(indPos(kk));
        dTerm(indPos(kk)) = -(sigmap*cosh(nuInfp*z)...
            +kappap*cosh(nuInfp*(z-1)))*csch(nuInfp)/...
            (sigmap*kappap*nuInfp);
      end
    end
    dTerm(isnan(dTerm)) = 0; % probably Rdl = Rf = 0, nonphysical short
  end
end