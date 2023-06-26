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
socfPct = 99;   % final cell SOC [%]
Crate = 1;      % discharge current [C-rate]
TdegC = 25;     % cell temperature [degC]
W = logspace(-1,1,10);  % Warburg resistance factor [unitless]
Wtau = 5;               % Warburg resistance factor when varying tau
taus = 0.1:0.2:0.7;     % Warburg time-constant for sep [s]
taup = 10:10:40;        % Warburg time-constant for pos [s]
Ts = min([taus taup])/10; % sampling interval [s]

% Load P2D cell model, convert to WORM.
p2dm = loadCellModel(modelName);
worm = convertCellModel(p2dm,'WORM');
Q = getCellParams(worm,'const.Q');

% Run simulations for different W.
clear WSeries;
for k = length(W):-1:1
    p.const.W = W(k);
    mod = setCellParam(worm,p);
    WSeries(k) = simCC( ...
        mod,Crate*Q,soc0Pct,socfPct,TdegC,'Ts',Ts,'Verbose',true);
end

% Run simulations for different taus.
clear tausSeries;
for k = length(taus):-1:1
    p.sep.tauW = taus(k);
    p.const.W = Wtau;
    mod = setCellParam(worm,p);
    tausSeries(k) = simCC( ...
        mod,Crate*Q,soc0Pct,socfPct,TdegC,'Ts',Ts,'Verbose',true);
end

% Run simulations for different taup.
clear taupSeries;
for k = length(taup):-1:1
    p.pos.tauW = taup(k);
    p.const.W = Wtau;
    mod = setCellParam(worm,p);
    taupSeries(k) = simCC( ...
        mod,Crate*Q,soc0Pct,socfPct,TdegC,'Ts',Ts,'Verbose',true);
end

simData.W = W;
simData.Wtau = Wtau;
simData.taus = taus;
simData.taup = taup;
simData.WSeries = WSeries;
simData.tausSeries = tausSeries;
simData.taupSeries = taupSeries;
simData.p2dm = p2dm;
simData.worm = worm;
simData.soc0Pct = soc0Pct;
simData.socfPct = socfPct;
simData.Crate = Crate;
simData.TdegC = TdegC;
save('warburgSim.mat','simData');
