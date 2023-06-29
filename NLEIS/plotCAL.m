% plotCAL.m
%
% Plot (NL)EIS calibration data!
%
% -- Changelog --
% 2023.05.11 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath(fullfile("..","UTILITY"));
expname = 'HybridNLEIS-CAL-Cell395534';
expdir = fullfile('labdata',expname);
plotdir = fullfile('plots',expname);
param = load(fullfile(expdir,'Parameters.mat'));
if ~isfolder(plotdir), mkdir(plotdir); end

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

% Plot support.
labels = arrayfun( ...
    @(x)sprintf('%.2f%s',x,param.cal_unit), ...
    param.cal_values,'UniformOutput',false);
labelNcol = 2;
colors = cool(length(param.cal_values));

% Nyquist plots.
figure;
colororder(colors);
plot(real(Z1),-imag(Z1),'.','MarkerFaceColor','auto');
legend(labels,'Location','best','NumColumns',labelNcol);
xlabel('Z_{1,1}'' [\Omega]');
ylabel('-Z_{1,1}'''' [\Omega]');
title('Nyquist: Linear Impedance');
setAxesNyquist;
thesisFormat([0.3 0.1 0.1 0.1]);
exportgraphics(gca,fullfile(plotdir,'Zcell-Nyq.eps'));
exportgraphics(gca,fullfile(plotdir,'Zcell-Nyq.png'));
figure;
colororder(colors);
plot(real(Z2),-imag(Z2),'.','MarkerFaceColor','auto');
xlabel('Z_{2,2}'' [\Omega A^{-1}]');
ylabel('-Z_{2,2}'''' [\Omega A^{-1}]');
title('Nyquist: Second-Harmonic Impedance');
legend(labels,'Location','best','NumColumns',labelNcol);
setAxesNyquist;
thesisFormat([0.3 0.1 0.1 0.1]);
exportgraphics(gca,fullfile(plotdir,'Zcell2-Nyq.eps'));
exportgraphics(gca,fullfile(plotdir,'Zcell2-Nyq.png'));

% Bode-magnitude plots.
figure;
colororder(colors);
loglog(freq,abs(Z1),'.','MarkerFaceColor','auto');
legend(labels,'Location','best','NumColumns',labelNcol);
xlabel('Cyclic Frequency [Hz]');
ylabel('|Z_{1,1}| [\Omega]');
title('Bode Magnitude: Linear Impedance');
thesisFormat([0.3 0.1 0.1 0.1]);
exportgraphics(gca,fullfile(plotdir,'Zcell-BodeMag.eps'));
exportgraphics(gca,fullfile(plotdir,'Zcell-BodeMag.png'));
figure;
colororder(colors);
loglog(freq,abs(Z2),'.','MarkerFaceColor','auto');
legend(labels,'Location','best','NumColumns',labelNcol);
xlabel('Cyclic Frequency [Hz]');
ylabel('|Z_{2,2}| [\Omega A^{-1}]');
title('Bode Magnitude: Second-Harmonic Impedance');
thesisFormat([0.3 0.1 0.1 0.1]);
exportgraphics(gca,fullfile(plotdir,'Zcell2-BodeMag.eps'));
exportgraphics(gca,fullfile(plotdir,'Zcell2-BodeMag.png'));

% Bode-phase plots.
figure;
colororder(colors);
semilogx(freq,(angle(Z1))*180/pi,'.','MarkerFaceColor','auto');
legend(labels,'Location','best','NumColumns',labelNcol);
xlabel('Cyclic Frequency [Hz]');
ylabel('\angle Z_{1,1} [deg]');
title('Bode Phase: Linear Impedance');
thesisFormat([0.3 0.1 0.1 0.1]);
exportgraphics(gca,fullfile(plotdir,'Zcell-BodePhs.eps'));
exportgraphics(gca,fullfile(plotdir,'Zcell-BodePhs.png'));
figure;
colororder(colors);
semilogx(freq,(angle(Z2))*180/pi,'.','MarkerFaceColor','auto');
legend(labels,'Location','best','NumColumns',labelNcol);
xlabel('Cyclic Frequency [Hz]');
ylabel('\angle Z_{2,2} [deg]');
title('Bode Phase: Second-Harmonic Impedance');
thesisFormat([0.3 0.1 0.1 0.1]);
exportgraphics(gca,fullfile(plotdir,'Zcell2-BodePhs.eps'));
exportgraphics(gca,fullfile(plotdir,'Zcell2-BodePhs.png'));

% KK Residuals.
figure;
colororder(colors);
semilogx(freq,100*abs(KK1resid));
legend(labels,'Location','best','NumColumns',labelNcol);
xlabel('Cyclic Frequency [Hz]');
ylabel('KK Residual, |Z_{kk}-Z||Z|^{-1} [%]');
title('KK Residual: Linear Impedance');
ylim([0 5]);
thesisFormat([0.3 0.1 0.1 0.1]);
exportgraphics(gca,fullfile(plotdir,'Zcell-KKresid.eps'));
exportgraphics(gca,fullfile(plotdir,'Zcell-KKresid.png'));

% Iapp THD.
figure;
colororder(colors);
semilogx(freq,Ithd);
legend(labels,'Location','best','NumColumns',labelNcol);
xlabel('Cyclic Frequency [Hz]');
ylabel('THD');
title('THD: Applied Current');
thesisFormat([0.3 0.1 0.1 0.1]);
exportgraphics(gca,fullfile(plotdir,'Iapp-THD.eps'));
exportgraphics(gca,fullfile(plotdir,'Iapp-THD.png'));

% V THD.
figure;
colororder(colors);
semilogx(freq,Vthd);
legend(labels,'Location','best','NumColumns',labelNcol);
xlabel('Cyclic Frequency [Hz]');
ylabel('THD');
title('THD: Cell Voltage');
thesisFormat([0.3 0.1 0.1 0.1]);
exportgraphics(gca,fullfile(plotdir,'Vcell-THD.eps'));
exportgraphics(gca,fullfile(plotdir,'Vcell-THD.png'));