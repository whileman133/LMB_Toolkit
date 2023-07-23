% plotFitEIS.m
%
% Plot regressed EIS model against lab data.

clear; close all; clc;
addpath(fullfile('..','..'));
TB.addpaths;
fitData = load(fullfile('labfitdata','EIS-Cell395524-42degC.mat'));
nsoc = length(fitData.socPctTrue);
indSOC = [1 2 3:2:nsoc-2 nsoc-1 nsoc];
nsocPlot = length(indSOC);
cellName = fitData.arg.labSpectra.cellName;
TdegC = fitData.TdegC;
plotdir = fullfile('plots',sprintf('LAB-Cell%s-%.0fdegC',cellName,TdegC));
if ~isfolder(plotdir)
    mkdir(plotdir);
end
rmseMag = sqrt(mean((abs(fitData.Zlab(:))-abs(fitData.Zmodel(:))).^2));
rmsePhs = sqrt(mean((angle(fitData.Zlab(:))-angle(fitData.Zmodel(:))).^2))*180/pi;
rmseReal = sqrt(mean((real(fitData.Zlab(:))-real(fitData.Zmodel(:))).^2));
rmseImag = sqrt(mean((imag(fitData.Zlab(:))-imag(fitData.Zmodel(:))).^2));

soc = linspace(0,1,100);
socPct = soc*100;
theta0 = fitData.values.pos.theta0;
theta100 = fitData.values.pos.theta100;
theta = theta0 + soc*(theta100-theta0);
electrode = MSMR(fitData.values.pos);
ctData = electrode.Rct(fitData.values.pos,'TdegC',TdegC,'theta',theta);
dsData = electrode.Ds(fitData.values.pos,'TdegC',TdegC,'theta',theta);

% Compare lab to fit impedance: Nyquist.
figure;
l = tiledlayout(3,ceil(nsocPlot/3));
l.Title.String = sprintf( ...
    'Linear Impedance Spectra: Cell %s (%.0f\\circC)',cellName,TdegC);
l.XLabel.String = 'Z'' [\Omega]';
l.YLabel.String = '-Z'''' [\Omega]';
for k = 1:nsocPlot
    ind = indSOC(k);
    zlab = fitData.Zlab(:,ind);
    zmod = fitData.Zmodel(:,ind);
    nexttile;
    plot(real(zlab),-imag(zlab),'b.'); hold on;
    plot(real(zmod),-imag(zmod),'r');
    title(sprintf('%.0f%% SOC',fitData.socPctTrue(ind)));
    setAxesNyquist;
end
legend('Lab','Fit Model','Location','northwest');
thesisFormat('FigSizeInches',[10 6]);
print(fullfile(plotdir,'Z-Nyq'),'-depsc');
print(fullfile(plotdir,'Z-Nyq'),'-dpng');

% Plot Rct.
figure;
semilogy(socPct,ctData.Rct,'k'); hold on;
semilogy( ...
    100*(fitData.values.pos.k0SplineTheta-theta0)/(theta100-theta0), ...
    1./fitData.values.pos.k0Spline./ctData.f,'ro');
ylabel('Charge-transfer resistance, R_{ct} [\Omega]');
xlabel('Cell state-of-charge [%]');
title(sprintf( ...
    'Log-Spline Charge-Transfer: Cell %s',cellName));
thesisFormat;
print(fullfile(plotdir,'Rct-soc'),'-depsc');
print(fullfile(plotdir,'Rct-soc'),'-dpng');

% Plot Ds.
figure;
semilogy(socPct,dsData.Ds,'k'); hold on;
semilogy( ...
    100*(fitData.values.pos.DsSplineTheta-theta0)/(theta100-theta0), ...
    fitData.values.pos.DsSpline,'ro');
ylabel('Solid diffusion coefficient, D_s [s^{-1}]');
xlabel('Cell state-of-charge [%]');
title(sprintf( ...
    'Log-Spline Solid Diffusion: Cell %s',cellName));
thesisFormat;
print(fullfile(plotdir,'Ds-soc'),'-depsc');
print(fullfile(plotdir,'Ds-soc'),'-dpng');