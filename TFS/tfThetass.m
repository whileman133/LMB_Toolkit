% function [thetassTF,aux] = tfThetass(s,locs,cellData) 
% 
% Inputs:
%   s         = vector of 1j*w where w is vector of frequencies at which
%               the transfer function (TF) is to be evaluated
%   locs      = vector of locations in normalized coordinates (0..3) at
%               which the TF is to be evaluated
%   cellData  = data structure containing all parameter values for a cell
%               at a given setpoint, most likely the output of
%               evalSetpoint.m 
% Outputs:
%   thetassTF = a matrix of the TF evaluated at every combination of "s"
%               and "locs" -- each row holds the TF for one location and
%               all frequencies; there is one row for every location
%   aux       = other outputs sometimes needed by calling routines (e.g.,
%               the xRAs). These outputs include dc-gains, high-frequency
%               gains, integrator residues, a record of the TF name, etc.
%
% This function evaluates the INTEGRATOR-REMOVED TF for normalized
% concentration of lithium in the solid, at the surface. If you wish to
% evaluate the TF with the integrator, use tfThetassInt.m instead.
%


function [thetassTF,aux] = tfThetass(s,locs,cellData) 
  % Initialize thetassstar to thetass
  s = s(:).';          % force "s" to be a row vector
  [C,Lambda,J,Z,Rct] = tfCommon(s,cellData); % Get Rct; store for efficiency
  cellData.common.C = C; cellData.common.L = Lambda; 
  cellData.common.J = J; cellData.common.Z = Z; cellData.common.Rct = Rct;
  [thetassTF,~] = tfThetassInt(s,locs,cellData);

  % Initialize outputs
  [aux.dcGain,r0p] = dcGain(locs);
  aux.hfGain = 0*locs(:);                % high-frequency gain  
  aux.res0 = 0*locs(:);                  % integrator residue
  aux.names = cell(length(locs),1);      % TF names
  aux.xLoc = zeros(1,length(locs));      % TF locations  


  indPos = find(locs >= 2);
  if ~isempty(indPos)
    aux.res0(indPos) = r0p;
    thetassTF(indPos,:) = thetassTF(indPos,:) ...
                              - r0p*ones(size(indPos(:)))*(1./s);
    [aux.names{indPos}] = deal('posThetass');
    aux.xLoc(indPos) = locs(indPos);                                                                                                    
  end
  inds0 = find(s==0,1);
  if ~isempty(inds0)
    thetassTF(:,inds0) = aux.dcGain;
  end

  function [dcGains,r0p] = dcGain(locs)
    dcGains = 0*locs(:); % dummy initialize
    
    T = cellData.const.T;
    
    % Set up constant values
    psi = cellData.const.psi;
    kappaD = cellData.const.kD;
    Q = cellData.const.Q; 
    


    % set up separator constants
    kappas = cellData.sep.kappa; % not SOC dependent
    qes = cellData.sep.qe;       % not SOC dependent

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
    Rctp = cellData.common.Rct(2);
    nDLp = cellData.pos.nDL;    
    Rdlp = cellData.pos.Rdl;
    wDLp = cellData.pos.wDL;
    CdlEffp = Cdlp^(2-nDLp)*wDLp^(nDLp-1);    
    r0p = -3*Dsp/(3*CdlEffp*Dsp*dUdqp-csmaxp); 

    s0 = 1e-7;
    [fr,~] = tfThetassInt(s0,locs,cellData);
    dcGains = fr - r0p*ones(size(fr(:)))*(1./s0);

    
      
    end    
  end
% end  