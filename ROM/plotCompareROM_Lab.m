% plotCompareROM_Lab.m

clear; close all; clc;
addpath('..');
TB.addpaths;

% Constants ---------------------------------------------------------------
simFile = 'Sion202309_GITT_Cell3955XX_P25_10m_60m_202310_EIS-16degC26degC-Ds=linear-k0=linear_defaultHRA';
data = load(fullfile('SIM_FILES',[simFile '.mat']));

ROMout = data.ROMout;
time = data.labData.time;
timeHr = time/60/60;
Vcell = data.labData.Vcell;
Iapp = data.labData.Iapp;

% Plotting ----------------------------------------------------------------
figure;
plot(timeHr,ROMout.Vcell); hold on;
plot(timeHr,Vcell,':');
xlim([min(timeHr) max(timeHr)]);
ylim([min(Vcell) max(Vcell)]);
title('Cell Voltage vs. Time: GITT (10m/60m 0.1C)');
xlabel('Time, t [hr]');
ylabel('v_{cell}(t)');
legend('HRA/Out-Blend Prediction','Laboratory Measurement','Location','northeast');
thesisFormat;
addInset([20 25],[350/60 3.3],2.5);

% RMSE Computation --------------------------------------------------------
ind = ~isnan(ROMout.Vcell(:));
vcellRMSE = sqrt(mean((ROMout.Vcell(ind)-Vcell(ind).').^2));
fprintf('RMSE:Vcell = %10.3f mV\n',vcellRMSE*1000);