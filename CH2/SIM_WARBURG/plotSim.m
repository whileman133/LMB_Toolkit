% plotSim.m
%
% Plot results of simulation study for Warburg resistance and time-constant 
% parameters.
%
% -- Changelog --
% 2023.06.23 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
load('warburgSim.mat');
labelsW = arrayfun( ...
    @(x)sprintf('$%.3f$',x),simData.W,'UniformOutput',false);
[Ld,Ls,Lp] = getCellParams(simData.p2dm,'*.L','Output','list');

figure; colororder(spring(length(simData.W)));
for k = 1:length(simData.W)
    data = simData.WSeries(k);
    Thetae = data.output.Thetae;  % dim1=time, dim2=xlocation
    xlocThetae = denormalizePosition(data.output.xLocs.Thetae,simData.p2dm);
    plot(1e6*xlocThetae,Thetae(end,:)); hold on;
end % for
xlim(1e6*[0 Ld+Ls+Lp]);
ylim(ylim());
rectangle('Position',[0, gca().YLim(1), Ld*1e6, diff(gca().YLim)],'FaceColor',[0 1 0 0.1],'LineStyle',':');
rectangle('Position',[Ld*1e6, gca().YLim(1), Ls*1e6, diff(gca().YLim)],'FaceColor',[0 0.5 1 0.1],'LineStyle',':');
rectangle('Position',[(Ld+Ls)*1e6, gca().YLim(1), Lp*1e6, diff(gca().YLim)],'FaceColor',[1 1 0 0.1],'LineStyle',':');
legend(labelsW{:},'Interpreter','latex','Location','best','NumColumns',5);
xlabel('Position along cell sandwich, $x$ [$\mu \mathrm{m}$]','Interpreter','latex');
ylabel('Salt concentration, $\theta_\mathrm{e}$ [unitless]','Interpreter','latex');
title('Pseudo-Steady Concentration Profiles vs. $\bar{W}$','Interpreter','latex');
thesisFormat;

figure; colororder(spring(length(simData.W)));
for k = 1:length(simData.W)
    data = simData.WSeries(k);
    Phie = data.output.Phie;  % dim1=time, dim2=xlocation
    xlocPhie = denormalizePosition(data.output.xLocs.Phie,simData.p2dm);
    plot(1e6*xlocPhie,Phie(end,:)); hold on;
end % for
xlim(1e6*[0 Ld+Ls+Lp]);
ylim(ylim());
rectangle('Position',[0, gca().YLim(1), Ld*1e6, diff(gca().YLim)],'FaceColor',[0 1 0 0.1],'LineStyle',':');
rectangle('Position',[Ld*1e6, gca().YLim(1), Ls*1e6, diff(gca().YLim)],'FaceColor',[0 0.5 1 0.1],'LineStyle',':');
rectangle('Position',[(Ld+Ls)*1e6, gca().YLim(1), Lp*1e6, diff(gca().YLim)],'FaceColor',[1 1 0 0.1],'LineStyle',':');
legend(labelsW{:},'Interpreter','latex','Location','best','NumColumns',5);
xlabel('Position along cell sandwich, $x$ [$\mu \mathrm{m}$]','Interpreter','latex');
ylabel('Electrolyte potential, $\phi_\mathrm{e}$ [V]','Interpreter','latex');
title('Pseudo-Steady Electrolyte Potential vs. $\bar{W}$','Interpreter','latex');
thesisFormat;

figure; colororder(spring(length(simData.taud)));
for k = 1:length(simData.taud)
    taud = simData.taud(k);
    data = simData.taudSeries(k);
    Phie = data.output.Phie;  % dim1=time, dim2=xlocation
    xlocPhie = data.output.xLocs.Phie;
    [~,indx] = min(abs(xlocPhie-1));
    indt = data.time<=max(simData.taud);
    plot(data.time(indt),Phie(indt,indx)); hold on;
end % for
thesisFormat;

figure; colororder(spring(length(simData.taus)));
for k = 1:length(simData.taus)
    taus = simData.taus(k);
    data = simData.tausSeries(k);
    Phie = data.output.Phie;  % dim1=time, dim2=xlocation
    xlocPhie = data.output.xLocs.Phie;
    [~,indx] = min(abs(xlocPhie-2));
    indt = data.time<=max(simData.taus);
    plot(data.time(indt),Phie(indt,indx)); hold on;
end % for
thesisFormat;

figure; colororder(spring(length(simData.taup)));
for k = 1:length(simData.taup)
    taup = simData.taup(k);
    data = simData.taupSeries(k);
    Phie = data.output.Phie;  % dim1=time, dim2=xlocation
    xlocPhie = data.output.xLocs.Phie;
    [~,indx] = min(abs(xlocPhie-3));
    indt = data.time<=max(simData.taup);
    plot(data.time(indt),Phie(indt,indx)); hold on;
end % for
thesisFormat;