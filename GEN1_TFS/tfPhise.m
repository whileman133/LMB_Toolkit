% function [phiseTF,aux] = tfPhise(s,locs,cellData) 
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
% This function evaluates the INTEGRATOR-REMOVED interphase-potential-
% difference (solid potential minus electrolyte potential, without
% integrator) TF. If you wish to evaluate the TF with the integrator, use
% tfPhiseInt.m instead.
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

function [phiseTF,aux] = tfPhise(s,locs,cellData)
  % Initialize phisestar to phise
  s = s(:).';          % force "s" to be a row vector
  [C,Lambda,J,Z,Rct] = tfCommon(s,cellData); % Get Rct; store for efficiency
  cellData.common.C = C; cellData.common.L = Lambda; 
  cellData.common.J = J; cellData.common.Z = Z; cellData.common.Rct = Rct;
  phiseTF = zeros(length(locs),length(s));
  aux.dcGain = zeros(length(locs),1);
  aux.hfGain = zeros(length(locs),1);

  aux.names = cell(length(locs),1);      % TF names
  aux.xLoc = zeros(1,length(locs));      % TF locations 
%   locs = locs';


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
  aux.res0 = 0*locs(indPos);                  % integrator residue

  if ~isempty(indPos)

   % Initialize outputs
    [phiseTFaux,aux2] = tfPhiseInt(s,locs(indPos),cellData);
    [aux.dcGain(indPos),r0p] = dcGain(locs(indPos));
 
    aux.hfGain(indPos) = aux2.hfGain;

    aux.res0(indPos) = r0p;
    phiseTF(indPos,:) = phiseTFaux ...
                            - r0p*ones(size(indPos(:)))*(1./s);
    [aux.names{indPos}] = deal('posPhise');
    aux.xLoc(indPos) = locs(indPos); 

     inds0 = find(s==0,1);
     if ~isempty(inds0)
         phiseTF(indPos,inds0) = aux.dcGain(indPos);  
     end
  end

  
  function [dcGains,r0p] = dcGain(locs)
    dcGains = 0*locs(:); % dummy initialize
    
    T = cellData.const.T;
    
    % Set up constant values
    
    kappaD = cellData.const.kD;
    Q = cellData.const.Q;
    
    % set up separator constants
    kappas = cellData.sep.kappa; 
    qes = cellData.sep.qe;       

    % set up positive-electrode constants...
    Dsp = cellData.pos.Ds;
    DeltaQp = abs(cellData.pos.theta100 - cellData.pos.theta0);
    csmaxp = 10800*Q*Dsp/DeltaQp;
    csmaxp2 = 3*Q*Dsp/DeltaQp;
    dUdqp = cellData.pos.dUocp;
    Cdlp = cellData.pos.Cdl;
    sigmap = cellData.pos.sigma;
    kappap = cellData.pos.kappa;
    qep = cellData.pos.qe;
    Rfp = cellData.pos.Rf; 
    Rctp = cellData.common.Rct(2);
    nDLp = cellData.pos.nDL;    
    Rdlp = cellData.pos.Rdl;
    wDLp = cellData.pos.wDL;
    CdlEffp = Cdlp^(2-nDLp)*wDLp^(nDLp-1);
    r0p = 3*Dsp*dUdqp/(csmaxp - 3*CdlEffp*Dsp*dUdqp);

    s0 = 1e-8;
    [fr,~] = tfPhiseInt(s0,locs,cellData);
    dcGains = fr - r0p*ones(size(fr(:)))*(1./s0);

    end
  end        
