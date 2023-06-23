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
% Values of warburg resistance factor for which to run the simulation.
W = logspace(-1,1,10);

% Load P2D cell model, convert to WORM.
p2dm = loadCellModel(modelName);
worm = convertCellModel(p2dm,'WORM');
Q = getCellParam(worm,'const.Q');

% Run simulations for different W.
clear WSeries;
for k = length(W):-1:1
    p.const.W = W(k);
    mod = setCellParam(worm,p);
    WSeries(k) = simCC(mod,Crate*Q,soc0Pct,socfPct,TdegC,'Verbose',true);
end

simData.W = W;
simData.WSeries = WSeries;
simData.p2dm = p2dm;
simData.worm = worm;
simData.soc0Pct = soc0Pct;
simData.socfPct = socfPct;
simData.Crate = Crate;
simData.TdegC = TdegC;
save('warburgSim.mat','simData');
