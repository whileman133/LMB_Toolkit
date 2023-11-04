% simFOM_MultipleHalfCycle.m
%
% Simulate slow half-cycle dis/charge response of the full-order LMB model 
% in COMSOL and save results to disk. Simulates over multiple C-rates.
%
% -- Changelog --
% 11.03.2023 | Incorporate chg response in addition to dischg | Wes H
% 09.07.2023 | Update for gen2 cell model | Wes H
% 05.31.2023 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath('..');
TB.addpaths;

% Constants.
cellFile = 'cellLMO-P2DM.xlsx';  % Name of cell parameters spreadsheet.
soc0Pct = 100;         % Initial cell SOC for dischg (final for chg) [%].
socfPct = 0;           % Final cell SOC for dischg (initial for chg) [%].
TdegC = 25;            % Cell temperature [degC].
IavgC = [-1:0.2:-0.6 0.6:0.2:1];   % Amplitude of Iapp sinusoid (- for chg) [average C-rate].
suffix = '';           % String to append to name of output data file.
% Structure of additional options to pass to simFOM.
OptSimFOM.VcellOnly = true;
OptSimFOM.DebugFlag = false;

% Load cell parameters.
cellModel = loadCellModel(cellFile);
cellModel = convertCellModel(cellModel,'RLWRM');
Q = getCellParams(cellModel,'const.Q');
Iavg = IavgC*Q; % average iapp in Amps

% Run EIS simulation(s) in COMSOL at each SOC setpoint.
clear iappSeries;
for k = length(Iavg):-1:1
    fprintf('Simulating Iavg=%5.2fC (#%d)... ',IavgC(k),k);
    zo = soc0Pct;
    zf = socfPct;
    if Iavg(k)<0
        % Discharge
        zo = socfPct;
        zf = soc0Pct;
    end
    iappSeries(k) = simHalfCycle( ...
        cellModel,abs(Iavg(k)),zo,zf,TdegC, ...
        'Verbose',false,'OptSimFOM',OptSimFOM);
    fprintf('done!\n')
end
simData.iappSeries = iappSeries;
simData.Iavg = Iavg;
simData.IavgC = IavgC;
simData.TdegC = TdegC;
simData.soc0Pct = soc0Pct;
simData.socfPct = socfPct;
simData.cellModel = cellModel;

% Save results to disk.
[~,modelName,~] = fileparts(cellFile);
fileName = fullfile( ...
    'simdata', ...
    'multhalfcyc', ...
    sprintf('%s-%dpct-%dpct',modelName,round(soc0Pct),round(socfPct)) ...
);
if ~isempty(suffix)
    fileName = [fileName '-' suffix];
end
save(fileName,"simData");