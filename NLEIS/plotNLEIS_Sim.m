% plotNLEIS_Sim.m

clear; close all; clc;
if ~exist('TOOLBOX_LMB','dir')
    % bootstrap the toolbox
    addpath('..');
    TB.addPaths();
end

% Process simulated spectra.
simName = 'cellLMO-Lumped-MSMR-30mA-socSeries';
load(fullfile('simdata',[simName '.mat']));
Z2sim = zeros(length(simData.freq),length(simData.socPct));
for k = 1:length(simData.socSeries)
    data = processEIS(simData.socSeries(k), 'NumHarmonics',2);
    Z2sim(:,k) = [data.h2.Zcell];
end

% Plotting ----------------------------------------------------------------
plotdir = fullfile('plots',['SIM_' simName]);
if ~isfolder(plotdir)
    mkdir(plotdir);
end
colors = cool(size(Z2sim,2));
labels = arrayfun( ...
    @(x)sprintf('%.0f%%',x),simData.socPct,'UniformOutput',false);

% Nyquist.
figure; colororder(colors);
plot(real(Z2sim),-imag(Z2sim),'.:','MarkerSize',15);
legend(labels{:},'NumColumns',2,'Location','best');
xlabel('$\tilde{Z}_{2,2}''\,[\Omega\,\mathrm{A}^{-1}]$','Interpreter','latex');
ylabel('-$\tilde{Z}_{2,2}''''\,[\Omega\,\mathrm{A}^{-1}]$','Interpreter','latex');
title(sprintf('Nyquist: $\\tilde{Z}_{2,2}$ (%s)', ...
    simData.cellModel.name),'Interpreter','latex');
setAxesNyquist('xdata',[-0.05, 0.05],'ydata',[-0.02, 0.02]);
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gca,fullfile(plotdir,'Z2-Nyq.eps'));
exportgraphics(gca,fullfile(plotdir,'Z2-Nyq.png'));