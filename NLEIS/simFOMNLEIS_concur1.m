% simFOMNLEIS.m
%
% Simulate medium-signal sinusoidal response of the full-order LMB model 
% in COMSOL and save results to disk.
%
% -- Changelog --
% 04.13.2023 | Allow using 'mesh' to specify all x-locs | Wesley Hileman
% 04.09.2023 | Develop simEIS.m utility function | Wesley Hileman
% 03.16.2023 | Fix issue in computing initial U0 | Wesley Hileman
% 03.09.2023 | Adapt for MSMR | Wesley Hileman
% 03.09.2023 | Increase sample rate 1 decade | Wesley Hileman
% 03.02.2023 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath(fullfile("..","UTILITY"));
addpath(fullfile("..","XLSX_CELLDEFS"));

% Constants.
cellFile = 'cellLMO-Lumped-fixDs.xlsx';  % Name of cell parameters spreadsheet.
freq = logspace(-3,5,50);   % Frequency points to evalulate in the spectrum [Hz].
socPct = 5;                 % Cell SOC setpoint [%].
TdegC = 25;                 % Cell temperature [degC].
I = 0.03;                   % Amplitude of Iapp sinusoids [A].
suffix = 'AllVars';         % String to append to name of output data file.
% Structure of additional options to pass to simFOM.
OptSimFOM.FixExchangeCurrent = true;

% The following structure specifies which electrochemical variables to
% store in addition to Vcell (field names) as well as the x-locations
% where those variables should be evaluated (field values).
% (Value of 'mesh' indicates all available x-locations should be stored.)
Vars.Phise = [0 3];
Vars.PhieTilde = 3;  % ground is at x=0+ (in the electrolyte)!
Vars.Thetae = 'mesh';
Vars.Thetass = 'mesh';
Vars.Eta = 'mesh';

% Load cell parameters and run EIS simulation in COMSOL.
cellModel = loadCellParams(cellFile);
simData = simEIS(cellModel,freq,socPct,TdegC, ...
    'Vars',Vars,'I',I,'Verbose',true,'OptSimFOM',OptSimFOM);

% Save results to disk.
[~,modelName,~] = fileparts(cellFile);
fileName = fullfile( ...
    'simdata', ...
    sprintf('%s-%dpct-%dmA', ...
        modelName,round(socPct),round(I*1000)) ...
);
if isfield(OptSimFOM,'FixExchangeCurrent') && OptSimFOM.FixExchangeCurrent
    fileName = [fileName '-fixI0'];
end
if ~isempty(suffix)
    fileName = [fileName '-' suffix];
end
save(fileName,"simData");