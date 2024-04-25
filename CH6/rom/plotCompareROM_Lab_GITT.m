% plotCompareROM_Lab_GITT.m

clear; close all; clc;
addpath('..');
TB.addpaths;

% Constants ---------------------------------------------------------------
simFileName = 'Sion202309_GITT_Cell3955XX_P25_10m_60m_EIS-16degC26degC-Ds=linear-k0=linear_defaultHRA';
data = load(fullfile('SIM_FILES',[simFileName '.mat']));

ROMout = data.ROMout;
time = data.labData.time;
timeHr = time/60/60;
Vcell = data.labData.Vcell;
Iapp = data.labData.Iapp;

% Plotting ----------------------------------------------------------------
plotdir = fullfile('PLOTS',simFileName);
if ~isfolder(plotdir)
    mkdir(plotdir);
end

figure;
plot(timeHr,ROMout.Vcell); hold on;
plot(timeHr,Vcell,':');
xlim([min(timeHr) max(timeHr)]);
ylim([min(Vcell) max(Vcell)]);
title('Vcell-v-Time');
xlabel('timehr');
ylabel('vcell');
legend('HRA/Out-Blend Prediction','Laboratory Measurement','Location','northeast');
thesisFormat;
addInset([20 25],[350/60 3.3],2.5);
print(fullfile(plotdir,'vcell'),'-depsc');
print(fullfile(plotdir,'vcell'),'-dpng');

% RMSE Computation --------------------------------------------------------
ind = ~isnan(ROMout.Vcell(:));
vcellRMSE = sqrt(mean((ROMout.Vcell(ind)-Vcell(ind).').^2));
fprintf('RMSE:Vcell = %10.3f mV\n',vcellRMSE*1000);