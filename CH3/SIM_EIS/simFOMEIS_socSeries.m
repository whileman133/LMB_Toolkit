% simFOMEIS_socSeries.m
%
% Simulate small-signal sinusoidal response of the full-order LMB model 
% in COMSOL and save results to disk. Simulates over multiple SOC
% setpoints.
%
% -- Changelog --
% 09.14.2023 | Update for gen2 toolkit | Wesley Hileman
% 05.30.2023 | Simulate over multiple SOC setpoints | Wesley Hileman
% 04.13.2023 | Allow using 'mesh' to specify all x-locs | Wesley Hileman
% 04.09.2023 | Develop simEIS.m utility function | Wesley Hileman
% 03.16.2023 | Fix issue in computing initial U0 | Wesley Hileman
% 03.09.2023 | Adapt for MSMR | Wesley Hileman
% 03.09.2023 | Increase sample rate 1 decade | Wesley Hileman
% 03.02.2023 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath(fullfile('..','..'));
TB.addpaths;

% Constants.
cellFile = 'cellLMO-P2DM.xlsx';  % Name of cell parameters spreadsheet.
freq = logspace(-4,5,50); % Frequency points to evalulate in the spectrum [Hz].
socPct = 10;               % Vector of cell SOC setpoint(s) [%].
TdegC = 25;               % Cell temperature [degC].
I = 0.001;                % Amplitude of Iapp sinusoids [A].
suffix = '';              % String to append to name of output data file.
% Structure of additional options to pass to simFOM.
OptSimFOM.FixExchangeCurrent = false;  % more realilistic when i0 varies
OptSimFOM.VcellOnly = false;

% The following structure specifies which electrochemical variables to
% store in addition to Vcell (field names) as well as the x-locations
% where those variables should be evaluated (field values).
% (Value of 'mesh' indicates all available x-locations should be stored.)
Vars.PhisTilde = 'mesh';  % debiased from vcell
Vars.Phise = 'mesh';
Vars.PhieTilde = 'mesh';  % ground is at x=0+ (in the electrolyte)
Vars.Eta = 'mesh';
Vars.Thetae = 'mesh';
Vars.Thetass = 'mesh';
Vars.Ifdl = 'mesh';
Vars.If = 'mesh';
Vars.Idl = 'mesh';

% Load cell parameters.
cellModel = loadCellModel(cellFile);

% Run EIS simulation(s) in COMSOL at each SOC setpoint.
clear socSeries;
for idxSOC = length(socPct):-1:1
    socSeries(idxSOC) = simEIS(cellModel,freq,socPct(idxSOC),TdegC, ...
        'Vars',Vars,'I',I,'Verbose',true,'OptSimFOM',OptSimFOM);
end
simData.socSeries = socSeries;
simData.socPct = socPct;
simData.TdegC = TdegC;
simData.I = I;
simData.cellModel = cellModel;

% Save results to disk.
[~,modelName,~] = fileparts(cellFile);
fileName = fullfile( ...
    'simdata', ...
    sprintf('%s-%dmA',modelName,round(I*1000)) ...
);
if isscalar(socPct)
    fileName = [fileName sprintf('-%dpct',round(socPct))];
else
    fileName = [fileName '-socSeries'];
end
if isfield(OptSimFOM,'FixExchangeCurrent') && OptSimFOM.FixExchangeCurrent
    fileName = [fileName '-fixI0'];
end
if ~isempty(suffix)
    fileName = [fileName '-' suffix];
end
save(fileName,"simData");