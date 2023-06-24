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
    @(x)sprintf('$\\bar{W}=%.3f$',x),simData.W,'UniformOutput',false);

figure; colororder(spring(length(simData.W)));
for k = 1:length(simData.W)
    W = simData.W(k);
    data = simData.WSeries(k);
    Thetae = data.output.Thetae;  % dim1=time, dim2=xlocation
    xlocThetae = denormalizePosition(data.output.xLocs.Thetae,simData.p2dm);
    plot(1e6*xlocThetae,Thetae(end,:)); hold on;
end % for
legend(labelsW{:},'Interpreter','latex','Location','best');
xlabel('Position along sandwich, $x$ [$\mu \mathrm{m}$]','Interpreter','latex');
ylabel('Salt concentration, $\theta_\mathrm{e}$ [unitless]','Interpreter','latex');
title('Pseudo-Steady-State Concentration Profiles');
thesisFormat;

figure; colororder(spring(length(simData.W)));
for k = 1:length(simData.W)
    W = simData.W(k);
    data = simData.WSeries(k);
    Phie = data.output.Phie;  % dim1=time, dim2=xlocation
    xlocPhie = denormalizePosition(data.output.xLocs.Phie,simData.p2dm);
    plot(1e6*xlocPhie,Phie(end,:)); hold on;
end % for
legend(labelsW{:},'Interpreter','latex','Location','best');
xlabel('Position along sandwich, $x$ [$\mu \mathrm{m}$]','Interpreter','latex');
ylabel('Electrolyte potential, $\phi_\mathrm{e}$ [V]','Interpreter','latex');
title('Pseudo-Steady-State Electrolyte Potential');
thesisFormat;