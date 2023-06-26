% runSim.m
%
% Run simulation study for layer reduction; combine dll and sep into single
% effective layer.
%
% -- Changelog --
% 2023.06.24 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath(fullfile('..','..'));
TB.addpaths('gen2');

modelName = 'cellLMO-P2DM';  % name of cell spreadsheet file
soc0Pct = 100;  % initial cell SOC [%]
socfPct = 98;   % final cell SOC [%]
Crate = 1;      % discharge current [C-rate]
TdegC = 25;     % cell temperature [degC]
Ts = 0.001;

% Load P2D cell model, convert to WORM.
p2dm = loadCellModel(modelName);
p.dll.eEps = 0.05;  % make dead-Li layer dense
p2dm = setCellParam(p2dm,p);
worm = convertCellModel(p2dm,'WORM');
rlworm = convertCellModel(worm,'RLWORM');
Q = getCellParams(worm,'const.Q');

% Run constant-current discharge simulations.
resultsWORM = simCC(worm,Crate*Q,soc0Pct,socfPct,TdegC,'Ts',Ts,'Verbose',true);
resultsRLWORM = simCC(rlworm,Crate*Q,soc0Pct,socfPct,TdegC,'Ts',Ts,'Verbose',true);

% Save results.
simData.resultsWORM = resultsWORM;
simData.resultsRLWORM = resultsRLWORM;
simData.p2dm = p2dm;
simData.worm = worm;
simData.rlworm = rlworm;
simData.soc0Pct = soc0Pct;
simData.socfPct = socfPct;
simData.Crate = Crate;
simData.TdegC = TdegC;
save('layerReductionSim.mat','simData');
