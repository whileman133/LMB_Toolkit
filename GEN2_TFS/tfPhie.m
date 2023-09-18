% function [phieTF,aux] = tfPhie(s,locs,cellData) 
% 
% Inputs:
%   s        = vector of 1j*w where w is vector of frequencies at which
%              the transfer function (TF) is to be evaluated
%   locs     = vector of locations in normalized coordinates (0..3) at
%              which the TF is to be evaluated
%   cellData = data structure containing all parameter values for a cell at
%              a given setpoint, most likely the output of evalSetpoint.m
% Outputs:
%   phieTF   = a matrix of the TF evaluated at every combination of "s" and
%              "locs" -- each row holds the TF for one location and all
%              frequencies; there is one row for every location
%   aux      = other outputs sometimes needed by calling routines (e.g.,
%              the xRAs). These outputs include dc-gains, high-frequency
%              gains, integrator residues, a record of the TF name, etc.
%
% This function evaluates the electrolyte-potential TF.
%
% Copyright (c) 2021 by Gregory L. Plett and M. Scott Trimboli of the
% University of Colorado Colorado Springs (UCCS). This work is licensed
% under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 Intl.
% License, v. 1.0. It is provided "as is", without express or implied
% warranty, for educational and informational purposes only.
%
% This file is provided as a supplement to: Plett, Gregory L. and Trimboli,
% M. Scott, "Battery Management Systems, Volume III, Physics-Based
% Methods," Artech House, 2021. 

function [phieTF,aux] = tfPhie(s,locs,cellData) 
  s = s(:).'; % force "s" to be a row vector
  [C,L,J,Z,Rct,param] = tfCommon(s,cellData); 
  cellData.common.C = C; 
  cellData.common.L = L; 
  cellData.common.J = J; 
  cellData.common.Z = Z;
  cellData.common.Rct = Rct;
  cellData.common.param = param;

  T = cellData.const.T;
  kappad = param.kappad;
  kappap = param.kappap;
  kappas = param.kappas;
  kappaD = param.kD;

  phieTF = zeros(length(locs),length(s));
  aux.hfGain = calcDterm(locs);          % high-frequency gain

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

  % Compute TF at each electrode
  indneg = find(locs == 0);
  if ~isempty(indneg)
      [fr, auxPhise] = tfPhise(s,0,cellData);
      phieTF(indneg,:) = -fr;
      aux.names{indneg} = 'negPhie';
      aux.xLoc(indneg) = 0;                          
      aux.hfGain(indneg) = -auxPhise.hfGain;
  end

  % [phien]1 & [phisen]2
  indDL = find(locs > 0 & locs <= 1);
  if ~isempty(indDL)
    for k = 1:length(indDL)
      z = locs(indDL(k));
      phieTF(indDL(k),:) = -z/kappad;
      phieTF(indDL(k),:) = phieTF(indDL(k),:) ...
          -kappaD*T*(C(c1DL,:).*(exp(L(L1DL,:)*(z-1))-exp(-L(L1DL,:))) ...
                   + C(c2DL,:).*(exp(-L(L1DL,:)*z)-1));
      aux.names{indDL(k)} = 'dlPhie';
      aux.xLoc(indDL(k)) = z;                          
    end
  end
  z = 1;
  phieDL1 = -z/kappad;
  phieDL1 = phieDL1  ...
          -kappaD*T*(C(c1DL,:).*(exp(L(L1DL,:)*(z-1))-exp(-L(L1DL,:))) ...
                   + C(c2DL,:).*(exp(-L(L1DL,:)*z)-1));

  % [phies]1 & [phies]2
  indSep = find(locs > 1 & locs < 2);
  if ~isempty(indSep)
    for k = 1:length(indSep)
      z = locs(indSep(k))-1;
      phieTF(indSep(k),:) = phieDL1  - z/kappas;
      phieTF(indSep(k),:) = phieTF(indSep(k),:) ...
          -kappaD*T*(C(c1s,:).*exp(L(L1s,:)*(z-1)) - C(c1s,:).*exp(-L(L1s,:)) ...
                   + C(c2s,:).*(exp(-L(L1s,:)*z)-1));
      aux.names{indSep(k)} = 'sepPhie';
      aux.xLoc(indSep(k)) = z+1;        
    end
  end
  z = 1;
  phies1 = phieDL1 - z/kappas;
  phies1 = phies1 ...
      -kappaD*T*(C(c1s,:).*exp(L(L1s,:)*(z-1)) - C(c1s,:).*exp(-L(L1s,:)) ...
               + C(c2s,:).*(exp(-L(L1s,:)*z)-1));

  % [phiep]1 & [phiep]2
  indPos = find(locs >= 2);
  if ~isempty(indPos)
    for k = 1:length(indPos)
      z = 3-locs(indPos(k));
      phieTF(indPos(k),:) = phies1  + (z-1)/kappap - ( ...
          +J(j1p,:).*(exp( L(L1p,:)*(z-1)) - 1           +(1-z)*L(L1p,:)                )./L(L1p,:).^2 ...
          +J(j2p,:).*(exp(-L(L1p,:)*z)    -exp(-L(L1p,:))-(1-z)*L(L1p,:).*exp(-L(L1p,:)))./L(L1p,:).^2 ...
          +J(j3p,:).*(exp( L(L2p,:)*(z-1))-1             +(1-z)*L(L2p,:)                )./L(L2p,:).^2 ...
          +J(j4p,:).*(exp(-L(L2p,:)*z)    -exp(-L(L2p,:))-(1-z)*L(L2p,:).*exp(-L(L2p,:)))./L(L2p,:).^2)/kappap;
      phieTF(indPos(k),:) = phieTF(indPos(k),:) ...
          - kappaD*T*(C(c1p,:).*(exp( L(L1p,:)*(z-1))-1             ) + ...
                     C(c2p,:).*(exp(-L(L1p,:)*z)    -exp(-L(L1p,:))) + ...
                     C(c3p,:).*(exp( L(L2p,:)*(z-1))-1             ) + ...
                     C(c4p,:).*(exp(-L(L2p,:)*z)    -exp(-L(L2p,:))));
      aux.names{indPos(k)} = 'posPhie';
      aux.xLoc(indPos(k)) = 3-z;                           
    end
  end


 % A function computes high-frequency gains (f->inf)
  function dTerm = calcDterm(locs)
    dTerm = 0*locs(:); % dummy initialize


    indDL = find(locs > 0 & locs <= 1);
    if ~isempty(indDL)
      for kk = 1:length(indDL)
        z = locs(indDL(kk));
        dTerm(indDL(kk)) = -z/kappad;
      end
    end

    indSep = find(locs > 1 & locs < 2);
    if ~isempty(indSep)
      for kk = 1:length(indSep)
        z = locs(indSep(kk))-1;
        dTerm(indSep(kk)) = -1/kappad -z/kappas;
      end
    end

    indPos = find(locs >= 2);
    Rf = cellData.pos.Rf;
    Rctp = cellData.common.Rct(2);
    ZseInfp = Rf + 1./(1./Rctp + 1./cellData.pos.Rdl);
    sigmap = cellData.pos.sigma;
    nuInfp = sqrt((1/sigmap+1/kappap)/ZseInfp);
    if ~isempty(indPos)
      for kk = 1:length(indPos)
        z = 1;
        dTerm(indPos(kk)) = -1/kappad -1/kappas;

        z = 3-locs(indPos(kk));
        dTerm(indPos(kk)) = dTerm(indPos(kk)) ...
            -(sigmap*(cosh(nuInfp)-cosh(nuInfp*z))...
            +kappap*(1-cosh(nuInfp*(1-z))+(1-z)*sinh(nuInfp)*nuInfp))...
            /(kappap*(sigmap+kappap)*sinh(nuInfp)*nuInfp);
      end
    end
  end
end