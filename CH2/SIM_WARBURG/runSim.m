% runSim.m
%
% Run simulation study for Warburg resistance and time-constant parameters.
%
% -- Changelog --
% 2023.06.23 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath(fullfile('..','..'));
TB.addpaths('gen2');

modelName = 'cellLMO-P2DM';  % name of cell spreadsheet file
soc0Pct = 100;  % initial cell SOC [%]
socfPct = 80;   % final cell SOC [%]
Crate = 1;      % discharge current [C-rate]
TdegC = 25;     % cell temperature [degC]

p2dm = loadCellModel(modelName);
worm = convertCellModel(p2dm,'WORM');
Q = getCellParam(worm,'const.Q');
simData = simCC(worm,Crate*Q,soc0Pct,socfPct,TdegC,'Verbose',true);
