% simGITTPulse.m
%
% Simulate relaxation of cell after application of a current pulse, both
% constant-current and quarter-sine, for single step of a GITT experiment.
%
% -- Changelog --
% 2023.03.10 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;

% Constants.
cellModelFile = 'cellLMO-P2DM.xlsx';  % cell model
soc0Pct = 100;    % initial SOC [%]
socfPct = 90;     % final SOC [%]
IavgCrate = 0.1;  % average dis/charge current [C-rate]
TdegC = 25;       % cell temperature [degC]
TrelaxSec = 2*3600; % amount of time to relax after pulse removed [s]

% Load cell model.
P2DM = loadCellModel(cellModelFile);
WRM = convertCellModel(P2DM,'WRM');
Q = getCellParams(WRM,'const.Q');
Iavg = IavgCrate*Q;

% Run simulations.
ccData = simCC(WRM,Iavg,soc0Pct,socfPct,TdegC,'Trelax',TrelaxSec);
qcData = simQuarterCycle(WRM,Iavg,soc0Pct,socfPct,TdegC,'Trelax',TrelaxSec);

% Save output.
save(fullfile('simdata','GITTPulse.mat'));