% demoSimEIS.m
% 
% Demonstrates the usage of the simEIS.m utility function for running
% full-order EIS simulations in COMSOL. Run this file in COMSOL with MATLAB
% LiveLink.
%
% Type `help simEIS` for detailed usage information.
%
% -- Changelog --
% 2023.04.05 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
thisfile = mfilename('fullpath');
[thisdir,thisfile,thisext] = fileparts(thisfile);
addpath(fullfile(thisdir,"TFS"));
addpath(fullfile(thisdir,"UTILITY"));
addpath(fullfile(thisdir,"XLSX_CELLDEFS"));

% Constants.
cellFile = 'cellLMO-Lumped-MSMR.xlsx';  % Name of cell parameters spreadsheet.
freq = logspace(-3,5,20);   % Frequency points to evalulate in the spectrum [Hz].
socPct = 50;                % Cell SOC setpoint [%].
TdegC = 25;                 % Cell temperature [degC].
I = 0.1;                    % Amplitude of Iapp sinusoids [A].
% Structure of additional options to pass to simFOM.
OptSimFOM.FixExchangeCurrent = true;

% The following structure specifies which electrochemical variables to
% store in addition to Vcell (field names) as well as the x-locations
% where those variables should be evaluated (field values).
Vars.Phise = [0 3];
Vars.Phie = 3;
Vars.Thetae = 0:3;

% Load cell parameters and run EIS simulation in COMSOL.
cellModel = loadCellParams(cellFile);
simData = simEIS(cellModel,freq,socPct,TdegC, ...
    'Vars',Vars,'I',I,'Verbose',true,'OptSimFOM',OptSimFOM);

% Save results to disk.
save("demoEIS.mat","simData");