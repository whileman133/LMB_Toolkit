% plotCompareROM_Lab_GITT.m

clear; close all; clc;
addpath('..');
TB.addpaths;

% Constants ---------------------------------------------------------------
romFilename = 'EIS-16degC26degC-Ds=linear-k0=linear_defaultHRA';
labDirectory = fullfile(TB.const.GITTROOT,'labdata');
labFilename = 'Sion202309_GITT_Cell3955XX_P25_10m_60m';
TdegC = 25;

% Collect lab data --------------------------------------------------------
file = fullfile(labDirectory,labFilename,'PWRGITT.DTA');
labData = loadGamryDTA(file,'NotesMode','KeyValue');
labTab = labData.tables.CURVE;
time = labTab.Time(:)';
Vcell = labTab.Vf(:)';
Iapp = -labTab.Im(:)'; % Gramy uses opposite sign convention
soc0Pct = labData.notes.soc0Pct;
name = labData.notes.name;

% Collect ROM -------------------------------------------------------------
romData = load(fullfile('ROM_FILES',[romFilename '.mat']));
ROM = romData.ROM;
LLPM = romData.LLPM;
[theta0,theta100] = getCellParams(ROM.cellData,'pos.theta0 pos.theta100','Output','list');
ocp0 = MSMR(ROM.cellData.function.pos).ocp('voltage',Vcell(1));
soc0PctTrue = 100*(ocp0.theta-theta0)/(theta100-theta0);

% Perform ROM simulation --------------------------------------------------
simData.SOC0 = soc0PctTrue;
simData.Iapp = Iapp;
simData.T = TdegC*ones(size(time));
simData.time = time;
ROMout = simROM(ROM,simData,'outBlend');

% Save results ------------------------------------------------------------
labData.time = time;
labData.Iapp = Iapp;
labData.Vcell = Vcell;
labData.soc0Pct = soc0Pct;
labData.soc0PctTrue = soc0PctTrue;
labData.name = name;
save(fullfile('SIM_FILES',[labFilename '_' romFilename '.mat']), ...
    'ROMout','romData','labData');