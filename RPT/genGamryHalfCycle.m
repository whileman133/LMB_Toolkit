% genGamryHalfCycle.m

clear; close all; clc;
addpath('..');
TB.addpaths;

Q = 0.06;
Iavg = Q/4;
soc0Pct = 100;
socfPct = 10;
Ts = 1;
name = 'HalfCyc_0C25_100pct_10pct_1s';
simData = simHalfCycle(Q,Iavg,soc0Pct,socfPct,'Ts',Ts,'DryRun',true);

figure;
plot(simData.time,simData.iapp);
thesisFormat;

tab = table(-simData.iapp);  % Gamry uses opposite sign convention!
writetable(tab,fullfile('gamry',[name '.txt']), ...
    'WriteVariableNames',false,'WriteRowNames',false);


