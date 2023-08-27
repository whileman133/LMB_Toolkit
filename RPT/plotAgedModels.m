% plotAgedModels.m
%
% Plot EIS, half-cycle discharge for various aged LMB models.
%
% -- Changelog --
% 2023.06.06 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath('..');
TB.addpaths();

cellDef = load('cellLMO_AgeSeries.mat');
TdegC = 25;
soc0Pct = 100;
socfPct = 0;
Iavg = 0.1;
freq = logspace(-3.5,6,100);

colors = spring(length(cellDef.lumped));
labels = arrayfun( ...
    @(x)sprintf('%.0f%%',x*100),cellDef.ageVect,'UniformOutput',false);

figure; colororder(colors);
for k = 1:length(cellDef.lumped)
    model = cellDef.lumped(k);
    halfcyc = simHalfCycle(model,Iavg,soc0Pct,socfPct,'DryRun',true);
    pertData = getPerturbationResistance(model,halfcyc.thetaAvg,'TdegC',TdegC);
    vcell = pertData.U - halfcyc.iapp.*pertData.Rtotal;
    plot(halfcyc.time/3600,vcell); hold on;
end % for
xlabel('Time, t [hr]');
ylabel('Terminal Voltage, v_{cell}(t) [V]');
title('Half-Cycle Discharge vs. State-of-Age');
legend(labels{:},'NumColumns',2);
thesisFormat([0.2 0.1 0.1 0.1]);

figure; colororder(colors);
for k = 1:length(cellDef.lumped)
    model = cellDef.lumped(k);
    tf = tfLMB(1j*2*pi*freq,model,'TdegC',TdegC,'socPct',5);
    Z1 = tf.h11.tfVcell();
    Z1 = Z1 - Z1(end);  % Subtract Z(inf)
    plot(real(Z1),-imag(Z1)); hold on;
end % for
xlabel('$[Z(f)-Z(\infty)]''$ [$\Omega$]','Interpreter','latex');
ylabel('$-[Z(f)-Z(\infty)]''''$ [$\Omega$]','Interpreter','latex');
title('Impedance Spectrum vs. State-of-Age (5% SOC)');
legend(labels{:},'NumColumns',2);
setAxesNyquist;
thesisFormat([0.2 0.1 0.1 0.1]);