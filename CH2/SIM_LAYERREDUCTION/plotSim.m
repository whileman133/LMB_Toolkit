% plotSim.m
%
% Plot results of simulation study for layer reduction.
%
% -- Changelog --
% 2023.06.24 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath(fullfile('..','..'));
TB.addpaths('gen2');
load('layerReductionSim.mat');

timeDesiredPlot = [0.1 62];  % times at which to plot conc. profiles
[Ld,Ls,Lp] = getCellParams(simData.p2dm,'*.L','Output','list');

time = simData.resultsWORM.time;
xtilde = simData.resultsWORM.output.xLocs.Phie;
xWORM = denormalizePosition(xtilde,simData.p2dm);
p.dll.L = (Ld+Ls)/2;
p.sep.L = (Ld+Ls)/2;
tmp = setCellParam(simData.p2dm,p);
xRLWORM = denormalizePosition(xtilde,tmp);
PhieWORM = simData.resultsWORM.output.Phie;
PhieRLWORM = simData.resultsRLWORM.output.Phie;
ThetaeWORM = simData.resultsWORM.output.Thetae;
ThetaeRLWORM = simData.resultsRLWORM.output.Thetae;

colors = summer(10);
colors = colors([1 6],:);
figure; colororder(colors);
plot(simData.resultsWORM.time, simData.resultsWORM.vcell,'-'); hold on;
plot(simData.resultsRLWORM.time, simData.resultsRLWORM.vcell,'--');
xlabel('Time, $t$ [s]','Interpreter','latex');
ylabel('Cell-level potential, $v_{\mathrm{cell}}$ [V]','Interpreter','latex');
title('Pulse-Response Simulation: Cell Voltage');
legend('dll and sep','eff only');
thesisFormat;
addInset([0 0.3],[10 3.237],2.7);
addInset([62 63],[48 3.22],2.7);
print('vcell','-depsc');
print('vcell','-dpng');

colors = spring(10);
colors = colors([2 9],:);
for t0 = timeDesiredPlot
    [~,indtime] = min(abs(time-t0));
    figure; colororder(colors);
    plot(xWORM*1e6,PhieWORM(indtime,:),'-'); hold on;
    plot(xRLWORM*1e6,PhieRLWORM(indtime,:),':');
    xlim([0 Ld+Ls+Lp]*1e6);
    ylim([-1.02 -0.94]);
    r1 = rectangle( ...
        'Position',[0, gca().YLim(1), Ld*1e6, diff(gca().YLim)], ...
        'FaceColor',[0 1 0 0.1], ...
        'LineStyle','-');
    r2 = rectangle( ...
        'Position',[Ld*1e6, gca().YLim(1), Ls*1e6, diff(gca().YLim)], ...
        'FaceColor',[0 0.5 1 0.1], ...
        'LineStyle','-');
    r3 = rectangle( ...
        'Position',[(Ld+Ls)*1e6, gca().YLim(1), Lp*1e6, diff(gca().YLim)], ...
        'FaceColor',[1 1 0 0.1], ...
        'LineStyle','-');
    uistack([r1 r2 r3],'bottom');
    set(gca, 'layer', 'top');
    xlabel('Position along cell sandwich, $x$ [$\mu \mathrm{m}$]','Interpreter','latex');
    ylabel('Electrolyte potential, $\phi_\mathrm{e}$ [V]','Interpreter','latex');
    title(sprintf('Electrolyte-Potential Profiles at $t=%.2f\\,\\mathrm{s}$',t0),'Interpreter','latex');
    legend('dll and sep','eff only');
    thesisFormat;
    print(sprintf('phie_%dms',t0*1000),'-depsc');
    print(sprintf('phie_%dms',t0*1000),'-dpng');
end
for t0 = timeDesiredPlot
    [~,indtime] = min(abs(time-t0));
    figure; colororder(colors);
    plot(xWORM*1e6,ThetaeWORM(indtime,:),'-'); hold on;
    plot(xRLWORM*1e6,ThetaeRLWORM(indtime,:),':');
    xlim([0 Ld+Ls+Lp]*1e6);
    ylim([0.95 1.1]);
    r1 = rectangle( ...
        'Position',[0, gca().YLim(1), Ld*1e6, diff(gca().YLim)], ...
        'FaceColor',[0 1 0 0.1], ...
        'LineStyle','-');
    r2 = rectangle( ...
        'Position',[Ld*1e6, gca().YLim(1), Ls*1e6, diff(gca().YLim)], ...
        'FaceColor',[0 0.5 1 0.1], ...
        'LineStyle','-');
    r3 = rectangle( ...
        'Position',[(Ld+Ls)*1e6, gca().YLim(1), Lp*1e6, diff(gca().YLim)], ...
        'FaceColor',[1 1 0 0.1], ...
        'LineStyle','-');
    uistack([r1 r2 r3],'bottom');
    set(gca, 'layer', 'top');
    xlabel('Position along cell sandwich, $x$ [$\mu \mathrm{m}$]','Interpreter','latex');
    ylabel('Salt concentration, $\theta_\mathrm{e}$ [unitless]','Interpreter','latex');
    title(sprintf('Concentration Profiles at $t=%.2f\\,\\mathrm{s}$',t0),'Interpreter','latex');
    legend('dll and sep','eff only');
    thesisFormat;
    print(sprintf('thetae_%dms',t0*1000),'-depsc');
    print(sprintf('thetae_%dms',t0*1000),'-dpng');
end