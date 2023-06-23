% function simData = loadInput(fileName)
% 
% Inputs:
%   fileName = filename of Excel spreadsheet to be loaded
%
% Outputs:
%   simData  = profile of current and temperature versus time; init SOC
%
% This function loads an Excel spreadsheet that specifies current and
% temperature versus time (plus initial SOC) for a ROM/FOM simulation.
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

function simData = loadInput(fileName)
  [num,~,~] = xlsread(fileName);
  simData.time = num(:,1);
  simData.time(isnan(num(:,1))) = [];
  simData.Iapp = num(:,2);
  simData.Iapp(isnan(num(:,2))) = [];
  simData.T    = num(:,3);
  simData.T(isnan(num(:,3))) = [];
  simData.Ts   = abs(num(2,1)-num(1,1));
  simData.SOC0 = num(:,4:5);
  simData.SOC0(isnan(num(:,4:5))) = [];
end
