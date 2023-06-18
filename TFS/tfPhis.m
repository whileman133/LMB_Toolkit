% function [phisTF,aux] = tfPhis(s,locs,cellData) 
% 
% Inputs:
%   s        = vector of 1j*w where w is vector of frequencies at which
%              the transfer function (TF) is to be evaluated
%   locs     = vector of locations in normalized coordinates (0..3) at
%              which the TF is to be evaluated
%   cellData = data structure containing all parameter values for a cell at
%              a given setpoint, most likely the output of evalSetpoint.m
% Outputs:
%   phisTF   = a matrix of the TF evaluated at every combination of "s" and
%              "locs" -- each row holds the TF for one location and all
%              frequencies; there is one row for every location
%   aux      = other outputs sometimes needed by calling routines (e.g.,
%              the xRAs). These outputs include dc-gains, high-frequency
%              gains, integrator residues, a record of the TF name, etc.
%
% This function evaluates the solid-potential TF.


function [phisTF,aux] = tfPhis(s,locs,cellData)
  % Set up some constants
  F = cellData.const.F; %#ok<*NASGU>
  T = cellData.const.T;

  % Initialize shared variables needed for calculating gains

  sigmap = cellData.pos.sigma;
  
  % Initialize outputs
  s = s(:).';          % force "s" to be a row vector
  [~,L,J,~,Rct] = tfCommon(s,cellData); 

  phisTF = zeros(length(locs),length(s));
  aux.hfGain = calcDterm(locs);          % high-frequency gain
  aux.dcGain = dcGain(locs);             % low-frequency gain
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



  indPos = find(locs >= 2);
  if ~isempty(indPos)
    for k = 1:length(indPos)
      z = 3-locs(indPos(k));
      phisTF(indPos(k),:) = (z ...
          +J(j1p,:).*(exp( L(L1p,:)*(z-1))-exp(-L(L1p,:)) ...
                  - z*L(L1p,:).*exp(-L(L1p,:)))./L(L1p,:).^2 ...
          +J(j2p,:).*(exp(-L(L1p,:)*z) - 1 ...
                  + z*L(L1p,:))./L(L1p,:).^2 ...
          +J(j3p,:).*(exp( L(L2p,:)*(z-1))-exp(-L(L2p,:)) ...
                  - z*L(L2p,:).*exp(-L(L2p,:)))./L(L2p,:).^2 ...
          +J(j4p,:).*(exp(-L(L2p,:)*z) - 1 ...
                  + z*L(L2p,:))./L(L2p,:).^2)/sigmap;
      aux.names{indPos(k)} = 'posPhis';
      aux.xLoc(indPos(k)) = 3-z;                                          
    end
  end

  if ~isempty(find(s==0,1))
      phisTF(:,s==0) = aux.dcGain;
  end

  % Computes DC gains (f->0)
  function dcGains = dcGain(locs)
    dcGains = 0*locs(:); % dummy initialize

    indPos = find(locs >= 2);
    if ~isempty(indPos)
      for kk = 1:length(indPos)
        z = 3-locs(indPos(kk));
        dcGains(indPos(kk)) = -z*(z-2)/(2*sigmap);
      end
    end
  end

  % Compute high-frequency gains (f->inf)
  function dTerm = calcDterm(locs)
    dTerm = 0*locs(:); % dummy initialize

    % set up positive-electrode constants...
    kappap = cellData.pos.kappa;
    Rfp = cellData.pos.Rf;
    Rctp = Rct(2);
    ZseInfp = Rfp + 1./(1./Rctp + 1./cellData.pos.Rdl);

    indPos = find(locs >= 2);
    if ~isempty(indPos)
      for kk = 1:length(indPos)
        z = 3-locs(indPos(kk));
        nuInfp = sqrt((1/sigmap+1/kappap)/ZseInfp);
        dTerm(indPos(kk)) = (kappap*(cosh(nuInfp) - cosh(nuInfp*(z-1))) +...
            sigmap*(1-cosh(z*nuInfp)+z*nuInfp*sinh(nuInfp)))...
            /(sigmap*(kappap+sigmap)*nuInfp*sinh(nuInfp));
      end
    end
  end
end