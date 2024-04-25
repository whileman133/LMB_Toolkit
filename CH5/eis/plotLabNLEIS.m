% plotLabNLEIS.m
%
% Plot (NL)EIS spectra collected in the lab!
%
% -- Changelog --
% 2023.05.15 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath('..');
TB.addpaths;
expname = '(NL)EIS-Cell395534_25degC';
expdir = fullfile(expname);
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
p = plot(real(spectra.lin.Z),-imag(spectra.lin.Z),'o-');
legend(labels{:},'NumColumns',4,'Location','southeast','FontSize',8);
xlabel('Re');
ylabel('-Im');
title('Zcell');
setAxesNyquist();
thesisFormat('LineMarkerSize',2,'LineMarkerLineWidth',0.5,'LineLineWidth',1);
ax = addInset([0.1 1.4],[5 5],'YSpan',[0 0.8]);
setAxesNyquist('axes',ax,'xdata',[0.1 1.4],'ydata',[0 0.8]);
print(fullfile(plotdir,'Zcell-Nyq'),'-depsc');
print(fullfile(plotdir,'Zcell-Nyq'),'-dpng');

% Bode: linear spectra.
figure; colororder(colors);
loglog(spectra.lin.freq,abs(spectra.lin.Z),'.');
xlabel('freq');
ylabel('MagZ');
title('BodeMag');
legend(labels{:},'NumColumns',3,'Location','best');
thesisFormat;
print(fullfile(plotdir,'Zcell-BodeMag'),'-depsc');
print(fullfile(plotdir,'Zcell-BodeMag'),'-dpng');
figure; colororder(colors);
semilogx(spectra.lin.freq,unwrap(angle(spectra.lin.Z))*180/pi,'.');
xlabel('freq');
ylabel('PhsZ');
title('BodePhs');
legend(labels{:},'NumColumns',3,'Location','best');
thesisFormat;
print(fullfile(plotdir,'Zcell-BodePhs'),'-depsc');
print(fullfile(plotdir,'Zcell-BodePhs'),'-dpng');