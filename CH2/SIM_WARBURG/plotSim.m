% plotSim.m
%
% Plot results of simulation study for Warburg resistance and time-constant 
% parameters.
%
% -- Changelog --
% 2023.06.23 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath(fullfile('..','..'));
TB.addpaths;

load('warburgSim.mat');
labelsW = arrayfun( ...
    @(x)sprintf('$\\bar{W}=%.3f$',x),simData.W,'UniformOutput',false);
labelsTaus = arrayfun( ...
    @(x)sprintf('$\\bar{\\tau}_\\mathrm{W}^\\mathrm{s}=%.3f$',x), ...
    simData.taus,'UniformOutput',false);
labelsTaup = arrayfun( ...
    @(x)sprintf('$\\bar{\\tau}_\\mathrm{W}^\\mathrm{p}=%.3f$',x), ...
    simData.taup,'UniformOutput',false);
[Ld,Ls,Lp] = getCellParams(simData.p2dm,'*.L','Output','list');

figure; colororder(copper(length(simData.W)));
for k = 1:length(simData.W)
    data = simData.WSeries(k);
    Thetae = data.output.Thetae;  % dim1=time, dim2=xlocation
    xlocThetae = denormalizePosition(data.output.xLocs.Thetae,simData.p2dm);
    plot(1e6*xlocThetae,Thetae(end,:)); hold on;
end % for
xlim(1e6*[0 Ld+Ls+Lp]);
ylim(ylim());
r1 = rectangle('Position',[0, gca().YLim(1), Ld*1e6, diff(gca().YLim)],'FaceColor',[0 1 0 0.1],'LineStyle','-');
r2 = rectangle('Position',[Ld*1e6, gca().YLim(1), Ls*1e6, diff(gca().YLim)],'FaceColor',[0 0.5 1 0.1],'LineStyle','-');
r3 = rectangle('Position',[(Ld+Ls)*1e6, gca().YLim(1), Lp*1e6, diff(gca().YLim)],'FaceColor',[1 1 0 0.1],'LineStyle','-');
uistack([r1 r2 r3],'bottom');
set(gca, 'layer', 'top');
legend(labelsW{:},'Interpreter','latex','Location','best','NumColumns',2);
xlabel('Position along cell sandwich, $x$ [$\mu \mathrm{m}$]','Interpreter','latex');
ylabel('Salt concentration, $\theta_\mathrm{e}$ [unitless]','Interpreter','latex');
title('Pseudo-Steady Concentration Profiles','Interpreter','latex');
thesisFormat;
print('thetaex-W','-depsc');
print('thetaex-W','-dpng');

figure; colororder(copper(length(simData.W)));
for k = 1:length(simData.W)
    data = simData.WSeries(k);
    Phie = data.output.Phie;  % dim1=time, dim2=xlocation
    xlocPhie = denormalizePosition(data.output.xLocs.Phie,simData.p2dm);
    plot(1e6*xlocPhie,Phie(end,:)); hold on;
end % for
xlim(1e6*[0 Ld+Ls+Lp]);
ylim(ylim());
r1 = rectangle('Position',[0, gca().YLim(1), Ld*1e6, diff(gca().YLim)],'FaceColor',[0 1 0 0.1],'LineStyle','-');
r2 = rectangle('Position',[Ld*1e6, gca().YLim(1), Ls*1e6, diff(gca().YLim)],'FaceColor',[0 0.5 1 0.1],'LineStyle','-');
r3 = rectangle('Position',[(Ld+Ls)*1e6, gca().YLim(1), Lp*1e6, diff(gca().YLim)],'FaceColor',[1 1 0 0.1],'LineStyle','-');
uistack([r1 r2 r3],'bottom');
set(gca, 'layer', 'top');
%legend(labelsW{:},'Interpreter','latex','Location','best','NumColumns',2);
xlabel('Position along cell sandwich, $x$ [$\mu \mathrm{m}$]','Interpreter','latex');
ylabel('Liquid-Phase potential, $\phi_\mathrm{e}$ [V]','Interpreter','latex');
title('Pseudo-Steady Electrolyte-Potential Profiles','Interpreter','latex');
thesisFormat;
print('phiex-W','-depsc');
print('phiex-W','-dpng');

figure; colororder(winter(length(simData.taus)));
for k = 1:length(simData.taus)
    taus = simData.taus(k);
    data = simData.tausSeries(k);
    Phie = data.output.Phie;  % dim1=time, dim2=xlocation
    xlocPhie = data.output.xLocs.Phie;
    [~,indx] = min(abs(xlocPhie-2));
    indt = data.time<=20*max(simData.taus);
    plot(data.time(indt),Phie(indt,indx)); hold on;
end % for
title('Potential at sep-pos Interface vs. Time', ...
    'Interpreter','latex');
xlabel('Time, $t$ [s]','Interpreter','latex');
ylabel('Liquid-Phase Potential, $\phi_\mathrm{e}$ [V]','Interpreter','latex');
legend(labelsTaus{:},'Interpreter','latex','Location','best');
thesisFormat;
print('phie2t-taus','-depsc');
print('phie2t-taus','-dpng');

figure; colororder(winter(length(simData.taup)));
for k = 1:length(simData.taup)
    taup = simData.taup(k);
    data = simData.taupSeries(k);
    Phie = data.output.Phie;  % dim1=time, dim2=xlocation
    xlocPhie = data.output.xLocs.Phie;
    [~,indx] = min(abs(xlocPhie-3));
    indt = data.time<=20*max(simData.taup);
    plot(data.time(indt),Phie(indt,indx)); hold on;
end % for
title('Potential at pos Current-Collector vs. Time',...
    'Interpreter','latex');
xlabel('Time, $t$ [s]','Interpreter','latex');
ylabel('Liquid-Phase Potential, $\phi_\mathrm{e}$ [V]','Interpreter','latex');
legend(labelsTaup{:},'Interpreter','latex','Location','best');
thesisFormat;
print('phie3t-taup','-depsc');
print('phie3t-taup','-dpng');