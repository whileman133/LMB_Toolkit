% plotLabNLEIS.m
%
% Plot (NL)EIS spectra collected in the lab!
%
% -- Changelog --
% 2023.05.15 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath(fullfile("..","UTILITY"));
expname = '(NL)EIS-SionCell395524';
expdir = fullfile('labdata',expname);
spectra = loadLabNLEIS(expdir);

colors = flip(cool(length(spectra.socPct)));
labels = arrayfun(@(x)sprintf('%d%%',x),spectra.socPct, ...
    'UniformOutput',false);
plotdir = fullfile('plots',expname);
if ~isfolder(plotdir)
    mkdir(plotdir);
end

% Nyquist: linear spectra.
figure; colororder(colors);
plot(real(spectra.lin.Z),-imag(spectra.lin.Z),'.:','MarkerSize',15);
legend(labels{:},'NumColumns',2,'Location','best');
xlabel('$\tilde{Z}_{1,1}''\,[\Omega]$','Interpreter','latex');
ylabel('-$\tilde{Z}_{1,1}''''\,[\Omega]$','Interpreter','latex');
title(sprintf('Nyquist: $\\tilde{Z}_{1,1}$ (Cell %s, %ddegC)',spectra.cellName,round(spectra.TdegC)), ...
    'Interpreter','latex');
setAxesNyquist( ...
    'xdata',[0 7], ...
    'ydata',[0 5]);
thesisFormat([0.2 0.1 0.2 0.1]);
exportgraphics(gcf,fullfile(plotdir,'Z1-Nyq.png'));
exportgraphics(gcf,fullfile(plotdir,'Z1-Nyq.eps'));

% Bode: linear spectra.
% figure; colororder(colors);
% loglog(spectra.lin.freq,abs(spectra.lin.Z),'.');
% legend(labels{:},'NumColumns',3,'Location','best');
% thesisFormat;
% figure; colororder(colors);
% semilogx(spectra.lin.freq,unwrap(angle(spectra.lin.Z))*180/pi,'.');
% legend(labels{:},'NumColumns',3,'Location','best');
% thesisFormat;

% Nyquist: nonlinear spectra.
idx = spectra.h2.freq<=10;
figure; colororder(colors);
plot(real(spectra.h2.Z(idx,:)),-imag(spectra.h2.Z(idx,:)),'.:','MarkerSize',15);
legend(labels{:},'NumColumns',3,'Location','best');
xlabel('$\tilde{Z}_{2,2}''\,[\Omega\,\mathrm{A}^{-1}]$','Interpreter','latex');
ylabel('-$\tilde{Z}_{2,2}''''\,[\Omega,\mathrm{A}^{-1}]$','Interpreter','latex');
title(sprintf('Nyquist: $\\tilde{Z}_{2,2}$ (Cell %s, %ddegC)',spectra.cellName,round(spectra.TdegC)), ...
    'Interpreter','latex');
setAxesNyquist;
thesisFormat([0.2 0.1 0.2 0.1]);
exportgraphics(gcf,fullfile(plotdir,'Z2-Nyq.png'));
exportgraphics(gcf,fullfile(plotdir,'Z2-Nyq.eps'));

% Bode: nonlinear spectra.
% figure; colororder(colors);
% loglog(spectra.h2.freq(idx),abs(spectra.h2.Z(idx,:)),'d:');
% legend(labels{:},'NumColumns',3,'Location','best');
% thesisFormat;
% figure; colororder(colors);
% semilogx(spectra.h2.freq(idx),unwrap(angle(spectra.h2.Z(idx,:)))*180/pi,'d:');
% legend(labels{:},'NumColumns',3,'Location','best');
% thesisFormat;