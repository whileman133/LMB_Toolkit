% function ROMout = simROM(ROM,simData,method,vmin,vmax,pmode)
% 
% Inputs:
%   ROM     = reduced-order model created with an XRA
%   simData = simulation profile loaded using "loadInput"
%   method  = switch between output/model/non blending ("outBlend",
%             "mdlBlend","nonBlend")
%   Vmin    = optional: minimum permitted cell voltage
%   Vmax    = optional: maximum permitted cell voltage
%   Pmode   = optional: if set to 1, input is power and not current
%
% Output:
%   ROMout  = data structure with results from simulation
%
% This utility function executes a reduced-order model using the
% user-prefered blending method
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

function ROMout = simROM(ROM,simData,method,varargin)
  Vmax = Inf; Vmin = -Inf; Pmode = 0;
  if nargin >= 4, Vmin  = varargin{1}; Vmax = varargin{2}; end
  if nargin >= 5, Pmode = varargin{3}; end

  % Update progress
  fprintf('----------------------------------------------------------\n');
  fprintf('Start to simulate ROM (%s)\n',datestr(now));

  % Simulate ROM using user defined blending methods
  if strcmpi(method,'outBlend')    
    ROMout = outBlend(simData,ROM,Vmin,Vmax,Pmode);
  elseif strcmpi(method,'mdlBlend')
    ROMout = mdlBlend(simData,ROM,Vmin,Vmax,Pmode);
  elseif strcmpi(method,'nonBlend')
    ROMout = nonBlend(simData,ROM,Vmin,Vmax,Pmode);
  elseif strcmpi(method,'nonBlendOld')
    ROMout = nonBlendOld(simData,ROM,Vmin,Vmax,Pmode);
  end

  fprintf('Finished: %s\n',datestr(now));
  fprintf('----------------------------------------------------------\n');
end