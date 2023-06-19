% processFOMNLEIS.m
%
% Plot results of medium-signal EIS simulation for full-order LMB cell.
%
% -- Changelog --
% 2023.04.05 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
thisfile = mfilename('fullpath');
[thisdir,thisfile,thisext] = fileparts(thisfile);
addpath(fullfile(thisdir,"TFS"));
addpath(fullfile(thisdir,"UTILITY"));
addpath(fullfile(thisdir,"XLSX_CELLDEFS"));
load('demoEIS.mat');

% Compute linear spectra from COMSOL data.
% Also evalulate linear TFs of same variables at same x-locs for 
% comparison to COMSOL simulation.
spectra = processEIS(simData, ...
    'NumHarmonics',1,'EvalLinTF',true,'NumTFFreqPoints',1000);
xxThetae = spectra.xlocs.Thetae;
xxPhise = spectra.xlocs.Phise;
xxPhie = spectra.xlocs.Phie;
Z1sim = [spectra.lin.Zcell];
Z1tf = [spectra.tf.Zcell];
ThetaeSim = [spectra.lin.Thetae];
ThetaeTF = [spectra.tf.Thetae];
PhiseSim = [spectra.lin.Phise];
PhiseTF = [spectra.tf.Phise];
PhieSim = [spectra.lin.Phie];
PhieTF = [spectra.tf.Phie];

% Plot linear impedance (Nyqiust).
figure;
plot(real(Z1tf),-imag(Z1tf)); hold on;
plot(real(Z1sim),-imag(Z1sim),'d');
title(sprintf('Nyquist: Z_{cell} (%s)',simData.param.cellModel.name));
xlabel('Z'' [\Omega]');
ylabel('-Z'''' [\Omega]');
legend('TF','COMSOL','Location','northwest');
setAxesNyquist;
thesisFormat;

% Plot linear impedance (Bode).
figure;
loglog(spectra.tfFreq,abs(Z1tf)); hold on;
loglog(spectra.freq,abs(Z1sim),'d');
xlabel('Cyclic Frequency [Hz]');
ylabel('|Z_{cell}| [\Omega]');
title(sprintf('Bode Magnitude: Z_{cell} (%s)',simData.param.cellModel.name));
legend('TF','COMSOL');
thesisFormat;
figure;
semilogx(spectra.tfFreq,angle(Z1tf)*180/pi); hold on;
semilogx(spectra.freq,angle(Z1sim)*180/pi,'d');
xlabel('Cyclic Frequency [Hz]');
ylabel('\angleZ_{cell} [deg]');
title(sprintf('Bode Phase: Z_{cell} (%s)',simData.param.cellModel.name));
legend('TF','COMSOL','Location','northwest');
thesisFormat;

% Plot linear Thetae/Iapp TF (Nyqiust).
labels1 = arrayfun(@(x)sprintf('x=%.0f TF',x),xxThetae, ...
    'UniformOutput',false);
labels2 = arrayfun(@(x)sprintf('x=%.0f FOM',x),xxThetae, ...
    'UniformOutput',false);
colors = cool(length(xxThetae));
figure;
for k = 1:length(xxThetae)
    plot(real(ThetaeTF(k,:)),-imag(ThetaeTF(k,:)), ...
        'Color',colors(k,:)); 
    hold on;
end
for k = 1:length(xxThetae)
    plot(real(ThetaeSim(k,:)),-imag(ThetaeSim(k,:)),'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
title('Nyquist: \Theta_e(j\omega)/Iapp(j\omega)');
xlabel('Real');
ylabel('-Imag');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','northwest');
setAxesNyquist;
thesisFormat([0.1 0.1 0.1 0.1]);

% Plot linear Phise/Iapp TF (Nyqiust).
labels1 = arrayfun(@(x)sprintf('x=%.0f TF',x),xxPhise, ...
    'UniformOutput',false);
labels2 = arrayfun(@(x)sprintf('x=%.0f FOM',x),xxPhise, ...
    'UniformOutput',false);
colors = cool(length(xxPhise));
figure;
for k = 1:length(xxPhise)
    plot(real(PhiseTF(k,:)),-imag(PhiseTF(k,:)), ...
        'Color',colors(k,:)); 
    hold on;
end
for k = 1:length(xxPhise)
    plot(real(PhiseSim(k,:)),-imag(PhiseSim(k,:)),'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
title('Nyquist: \Phi_{s,e}(j\omega)/Iapp(j\omega)');
xlabel('Real');
ylabel('-Imag');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','northwest');
setAxesNyquist;
thesisFormat([0.1 0.1 0.1 0.1]);

% Plot linear Phie/Iapp TF (Nyqiust).
labels1 = arrayfun(@(x)sprintf('x=%.0f TF',x),xxPhie, ...
    'UniformOutput',false);
labels2 = arrayfun(@(x)sprintf('x=%.0f FOM',x),xxPhie, ...
    'UniformOutput',false);
colors = cool(length(xxPhise));
figure;
for k = 1:length(xxPhie)
    plot(real(PhieTF(k,:)),-imag(PhieTF(k,:)), ...
        'Color',colors(k,:)); 
    hold on;
end
for k = 1:length(xxPhie)
    plot(real(PhieSim(k,:)),-imag(PhieSim(k,:)),'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
title('Nyquist: \Phi_{e}(j\omega)/Iapp(j\omega)');
xlabel('Real');
ylabel('-Imag');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','northwest');
setAxesNyquist;
thesisFormat([0.1 0.1 0.1 0.1]);