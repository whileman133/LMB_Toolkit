% function [idlTF,aux] = tfIdl(s,locs,cellData) 
% 
% Inputs:
%   s        = vector of 1j*w where w is vector of frequencies at which
%              the transfer function (TF) is to be evaluated
%   locs     = vector of locations in normalized coordinates (0..3) at
%              which the TF is to be evaluated
%   cellData = data structure containing all parameter values for a cell at
%              a given setpoint, most likely the output of evalSetpoint.m
% Outputs:
%   idlTF    = a matrix of the TF evaluated at every combination of "s" and
%              "locs" -- each row holds the TF for one location and all
%              frequencies; there is one row for every location
%   aux      = other outputs sometimes needed by calling routines (e.g.,
%              the xRAs). These outputs include dc-gains, high-frequency
%              gains, integrator residues, a record of the TF name, etc.
%
% This function evaluates the double-layer (nonfaradaic) lithium flux TF. 
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

function [idlTF,aux] = tfIdl(s,locs,cellData) %#ok<*NASGU>
  % Set up some constants
  F = cellData.const.F;
  Q = cellData.const.Q;

  % Initialize outputs
  s = s(:).';  % force "s" to be a row vector
  [C,L,J,Z,Rct] = tfCommon(s,cellData); % Get Rct; store for efficiency
  cellData.common.C = C; cellData.common.L = L; 
  cellData.common.J = J; cellData.common.Z = Z; cellData.common.Rct = Rct;
  
  idlTF = zeros(length(locs),length(s));
  aux.hfGain = calcDterm(locs);      % high-frequency gain
  aux.dcGain = dcGain(locs);         % low-frequency gain
  aux.res0 = 0*locs(:);              % integrator residue
  aux.names = cell(length(locs),1);  % TF names
  aux.xLoc = zeros(1,length(locs));  % TF locations

  % Get common components... (next three lines needed to keep MATLAB happy)
  L1DL=[];L1s=[];L1p=[];L2p=[];
  c1DL=[];c2DL=[];
  c1s=[];c2s=[];
  c1p=[];c2p=[];c3p=[];c4p=[];
  j1p=[];j2p=[];j3p=[];j4p=[];
  Zsep=[]; Zsp=[];
  Isn=[];Isp=[]; 
  eval(cellData.common.ind); % set c##,L##,J##,Z## indices


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
    Ze = 0;  
    Zf = Zdln.*(Rctn + Ze)./(Zdln - Ze);
%     Zf = Rctn;
    Zf = Zdln.*(Rctn + Ze)./(Zdln - Ze);
    idl = Zf./(Zdln + Zf);
    idlTF(indNeg,:) = idl(2:end);

    s = s(2:end);

    [aux.names{indNeg}] = deal('negIdl');
    aux.xLoc(indNeg) = locs(indNeg);  
    aux.hfGain(indNeg) = Rctn/(Rdln + Rctn);

     inds0 = find(s==0,1);
     if ~isempty(inds0)
        aux.dcGain(indNeg) = idlTF(indNeg,1);
        idlTF(indNeg,inds0) = aux.dcGain(indNeg); 

     end
  end



  % Compute Rctn and Rctp
  Rctp = cellData.common.Rct(2);
  Rdlp = cellData.pos.Rdl;
  Cdlp = cellData.pos.Cdl;
  nDLp = cellData.pos.nDL;
  wDLp = cellData.pos.wDL;  
  Zdlp = Rdlp + (Cdlp*s+wDLp).^(1-nDLp)./(Cdlp^(2-nDLp)*s);


  % Compute TF in positive electrode
  indPos = find(locs >= 2);
  if ~isempty(indPos)
    for k = 1:length(indPos)
      z = 3-locs(indPos(k));
      idl = ...
          J(j1p,:).*exp(L(L1p,:)*(z-1)) + J(j2p,:).*exp(-L(L1p,:)*z) + ...
          J(j3p,:).*exp(L(L2p,:)*(z-1)) + J(j4p,:).*exp(-L(L2p,:)*z);
      idlTF(indPos(k),:) = idl.*Z(Isp,:)./Zdlp.*(Rctp+Z(Zsp,:));
      aux.names{indPos(k)} = 'posIdl';
      aux.xLoc(indPos(k)) = 3-z;
    end
  end

  % Force dcGain
  if ~isempty(find(s==0,1))
%     idlTF(:,s==0) = aux.dcGain;
  end

  % Computes high-frequency gains (f->inf)
  function dTerm = calcDterm(locs)
    dTerm = 0*locs(:); % dummy initialize



    % set up positive-electrode constants...
    kappap = cellData.pos.kappa;
    Rfp = cellData.pos.Rf;
    Rctp = cellData.common.Rct(2);
    ZseInfp = Rfp + 1./(1./Rctp + 1./cellData.pos.Rdl);
    sigmap = cellData.pos.sigma;


    indPos = find(locs >= 2);
    if ~isempty(indPos)
      for kk = 1:length(indPos)
        z = 3-locs(indPos(kk));
        nuInfp = sqrt((1/sigmap+1/kappap)/ZseInfp);
        dTerm(indPos(kk)) = -nuInfp*(sigmap*cosh(nuInfp*z)+...
            kappap*cosh(nuInfp*(z-1)))...
            /((kappap+sigmap)*sinh(nuInfp))...
            *1/(cellData.pos.Rdl/Rctp+1);
      end
    end
  end

  % Computes dc gains (f->0)
  function dcGains = dcGain(locs)
    dcGains = 0*locs(:); % dummy initialize

    indNeg = find(locs == 0);
    dcGains(indNeg) = 0;

    % set up positive-electrode constants...
    dUdqp = cellData.pos.dUocp;
    theta0p = cellData.pos.theta0;
    theta100p = cellData.pos.theta100;
    thetapABS = abs(theta100p-theta0p);
    Cdlp = cellData.pos.Cdl;
    nDLp = cellData.pos.nDL;
    wDLp = cellData.pos.wDL;    
    CdlEffp = Cdlp^(2-nDLp)*wDLp^(nDLp-1);


    indPos = find(locs >= 2);
    if ~isempty(indPos)
      for kk = 1:length(indPos)
        z = 3-locs(indPos(kk));
        dcGains(indPos(kk)) = CdlEffp*thetapABS*dUdqp...
            /(3600*Q-CdlEffp*thetapABS*dUdqp);
      end
    end
  end
end