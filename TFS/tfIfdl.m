% function [ifdlTF,aux] = tfIfdl(s,locs,cellData) 
% 
% Inputs:
%   s        = vector of 1j*w where w is vector of frequencies at which
%              the transfer function (TF) is to be evaluated
%   locs     = vector of locations in normalized coordinates (0..3) at
%              which the TF is to be evaluated
%   cellData = data structure containing all parameter values for a cell at
%              a given setpoint, most likely the output of evalSetpoint.m
% Outputs:
%   ifdlTF   = a matrix of the TF evaluated at every combination of "s" and
%              "locs" -- each row holds the TF for one location and all
%              frequencies; there is one row for every location
%   aux      = other outputs sometimes needed by calling routines (e.g.,
%              the xRAs). These outputs include dc-gains, high-frequency
%              gains, integrator residues, a record of the TF name, etc.
%
% This function evaluates the total (faradaic plus nonfaradaic) lithium
% flux TF. 
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

function [ifdlTF,aux] = tfIfdl(s,locs,cellData) %#ok<*NASGU>
  % Initialize outputs
  s = s(:).';          % force "s" to be a row vector  
  [C,L,J,Z,Rct] = tfCommon(s,cellData); % Get Rct; store for efficiency
  cellData.common.C = C; cellData.common.L = L; 
  cellData.common.J = J; cellData.common.Z = Z; cellData.common.Rct = Rct;
  
  ifdlTF = zeros(length(locs),length(s));
  aux.hfGain = calcDterm(locs);          % high-frequency gain
  aux.dcGain = ones(length(locs),1);     % low-frequency gain
  aux.dcGain(locs>1.5) = -aux.dcGain(locs>1.5);
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
    [If0TF,auxIf] = tfIf(s,0,cellData);
    [Idl0TF,auxIdl] = tfIdl(s,0,cellData);

    ifdlTF(indNeg,:) = If0TF(2:end) + Idl0TF(2:end);


    s = s(2:end);

    [aux.names{indNeg}] = deal('negIfdl');
    aux.xLoc(indNeg) = locs(indNeg);  
    aux.hfGain(indNeg) = 1;

     inds0 = find(s==0,1);
     if ~isempty(inds0)
        aux.dcGain(indNeg) = If0TF(1) + Idl0TF(1);
        ifdlTF(indNeg,inds0) = aux.dcGain(indNeg); 

     end
  end



  indPos = find(locs >= 2);
  if ~isempty(indPos)
    for k = 1:length(indPos)
      z = 3-locs(indPos(k));
      ifdlTF(indPos(k),:) = ...
          J(j1p,:).*exp(L(L1p,:)*(z-1)) + J(j2p,:).*exp(-L(L1p,:)*z) + ...
          J(j3p,:).*exp(L(L2p,:)*(z-1)) + J(j4p,:).*exp(-L(L2p,:)*z);
      aux.names{indPos(k)} = 'posIfdl';
      aux.xLoc(indPos(k)) = 3-z;                
    end
  end

  if ~isempty(find(s==0,1))
    ifdlTF(:,s==0) = aux.dcGain;
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
            kappap*cosh(nuInfp*(z-1)))/((kappap+sigmap)*sinh(nuInfp));
      end
    end
  end
end