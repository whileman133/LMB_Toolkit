% function [thetassTF,aux] = tfThetassInt(s,locs,cellData) 
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
% This function evaluates the TF for normalized concentration of lithium in
% the solid, at the surface WITH the integration dynamics included. If you
% wish to evaluate the TF without the integrator, use tfThetass.m instead.
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

function [thetassTF,aux] = tfThetassInt(s,locs,cellData) %#ok<*NASGU>
  % Set up some constants
  F = cellData.const.F; 
  T = cellData.const.T;

  % Initialize outputs
  s = s(:).';          % force "s" to be a row vector
  thetassTF = zeros(length(locs),length(s));
%   aux.hfGain = 0*locs(:);                % high-frequency gain
%   aux.dcGain = Inf(length(locs),1);      % low-frequency gain
%   aux.res0 = 0*locs(:);                  % integrator residue
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

  [~,L,J,Z] = tfCommon(s,cellData); 
%   eval(cellData.common.ind); % set c##,L##,J##,Z## indices


  % Compute TF in positive electrode
  indPos = find(locs >= 2);  
  if ~isempty(indPos)
    dUdqp = cellData.pos.dUocp;    
    for k = 1:length(indPos)
      z = 3-locs(indPos(k));
      ifdl_tf = ...
          J(j1p,:).*exp(L(L1p,:)*(z-1)) + J(j2p,:).*exp(-L(L1p,:)*z) + ...
          J(j3p,:).*exp(L(L2p,:)*(z-1)) + J(j4p,:).*exp(-L(L2p,:)*z);
      thetassTF(indPos(k),:) = (1/dUdqp)*ifdl_tf.*Z(Isp,:).*Z(Zsp,:);
      aux.names{indPos(k)} = 'posThetassInt';
      aux.xLoc(indPos(k)) = 3-z;
    end
  end
end