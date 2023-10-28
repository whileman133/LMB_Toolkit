% plotCompareROM_Lab.m

clear; close all; clc;
addpath('..');
TB.addpaths;

% Constants ---------------------------------------------------------------
romFilename = '202309_EIS-16degC26degC-Ds=linear-k0=linear_defaultHRA';
labDirectory = fullfile(TB.const.GITTROOT,'labdata');
labFilename = 'Sion202309_GITT_Cell3955XX_P25_10m_60m';
TdegC = 25;

% Collect lab data --------------------------------------------------------
file = fullfile(labDirectory,labFilename,'PWRGITT.DTA');
labData = loadGamryDTA(file,'NotesMode','KeyValue');
labTab = labData.tables.CURVE;
time = labTab.Time(:)';
timeMin = time/60;
Vcell = labTab.Vf(:)';
Iapp = -labTab.Im(:)'; % Gramy uses opposite sign convention
soc0Pct = labData.notes.soc0Pct;
name = labData.notes.name;

% Collect ROM -------------------------------------------------------------
romData = load(fullfile('ROM_FILES',[romFilename '.mat']));
ROM = romData.ROM;
LLPM = romData.LLPM;

% Perform ROM simulation --------------------------------------------------
simData.SOC0 = soc0Pct;
simData.Iapp = Iapp;
simData.T = TdegC*ones(size(time));
simData.time = time;
ROMout = simROM(ROM,simData,'outBlend');

% Plotting ----------------------------------------------------------------
figure;
plot(timeMin,ROMout.Vcell); hold on;
plot(timeMin,Vcell,':');
title('Cell Voltage vs. Time');
xlabel('Time, t [min]');
ylabel('v_{cell}(t)');
legend('ROM','Lab','Location','best');
thesisFormat;

% RMSE Computation --------------------------------------------------------
vcellRMSE = sqrt(mean((ROMout.Vcell(:)-Vcell(:)).^2));
fprintf('RMSE:Vcell = %10.3f mV\n',vcellRMSE*1000);