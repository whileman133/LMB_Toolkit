% plotCAL.m
%
% Plot (NL)EIS calibration data!
%
% -- Changelog --
% 2023.05.11 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath(fullfile('..','..'));
TB.addpaths;

% Constanta.
expname = 'GalvNLEIS-CAL-Cell395534';
expdir = fullfile(expname);
plotdir = fullfile('plots',expname);

% Load experiment parameters.
param = load(fullfile(expdir,'Parameters.mat'));

% Collect impedance spectra.
Z1 = [];
Z2 = [];
KK1resid = [];
Ithd = [];
Vthd = [];
for k = size(param.eis_filenames,1):-1:1
    filename = strtrim(param.eis_filenames(k,:));
    eis = loadGamryEIS(fullfile(expdir,filename));
    KK1 = linKK(eis.freq,eis.Z,'IncludeIntegrator',true,'c',0.9);
    Z1(:,k) = eis.Zhh(:,1);
    Z2(:,k) = eis.Zhh(:,2);
    KK1resid(:,k) = KK1.residuals;
    Ithd(:,k) = eis.Ithd;
    Vthd(:,k) = eis.Vthd;
end
freq = eis.freq;

% Plotting ----------------------------------------------------------------
if ~isfolder(plotdir), mkdir(plotdir); end
labels = arrayfun( ...
    @(x)sprintf('%.2f%s',x,param.cal_unit), ...
    param.cal_values,'UniformOutput',false);
labelNcol = 2;
colors = cool(length(param.cal_values));

% Nyquist plots.
figure;
colororder(colors);
plot(real(Z1),-imag(Z1),'.','MarkerFaceColor','auto');
legend(labels,'Location','south','NumColumns',labelNcol);
xlabel('Re');
ylabel('-Im');
title('Zcell');
setAxesNyquist;
thesisFormat();
print(fullfile(plotdir,'Zcell-Nyq'),'-depsc');
print(fullfile(plotdir,'Zcell-Nyq'),'-dpng');

% KK Residuals.
figure;
colororder(colors);
semilogx(freq,100*abs(KK1resid));
xlim([min(freq) max(freq)]);
legend(labels,'Location','best','NumColumns',labelNcol);
xlabel('freq');
ylabel('PctResidual');
title('KKResid');
ylim([0 5]);
thesisFormat;
print(fullfile(plotdir,'Zcell-KKresid'),'-depsc');
print(fullfile(plotdir,'Zcell-KKresid'),'-dpng');

% Iapp THD.
figure;
colororder(colors);
semilogx(freq,Ithd*100);
xlim([min(freq) max(freq)]);
legend(labels,'Location','best','NumColumns',labelNcol);
xlabel('freq');
ylabel('THD');
title('THD-iapp');
thesisFormat;
print(fullfile(plotdir,'Iapp-THD'),'-depsc');
print(fullfile(plotdir,'Iapp-THD'),'-dpng');

% V THD.
figure;
colororder(colors);
semilogx(freq,Vthd*100);
xlim([min(freq) max(freq)]);
legend(labels,'Location','best','NumColumns',labelNcol);
xlabel('freq');
ylabel('THD');
title('THD-vcell');
thesisFormat;
print(fullfile(plotdir,'Vcell-THD'),'-depsc');
print(fullfile(plotdir,'Vcell-THD'),'-depsc');