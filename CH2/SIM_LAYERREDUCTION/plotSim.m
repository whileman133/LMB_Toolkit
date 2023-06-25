% plotSim.m
%
% Plot results of simulation study for layer reduction.
%
% -- Changelog --
% 2023.06.24 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
load('layerReductionSim.mat');

colors = spring(10);
colors = colors([2 9],:);
figure; colororder(colors);
plot(simData.resultsWORM.time, simData.resultsWORM.vcell,'-'); hold on;
plot(simData.resultsRLWORM.time, simData.resultsRLWORM.vcell,'--');
xlabel('Time, $t$ [s]','Interpreter','latex');
ylabel('Cell-level potential, $v_{\mathrm{cell}}$ [V]','Interpreter','latex');
title('Pulse-Response Simulation: Cell Voltage');
legend('dll and sep','eff only');
thesisFormat;
addInsetAxes([0 0.3],[10 3.266],2.7);
addInsetAxes([60 65],[48 3.248],2.7);
print('vcell','-depsc');
print('vcell','-dpng');