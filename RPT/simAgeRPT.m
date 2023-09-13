% simAgeRPT.m
%
% Simulate full RPT (half-cycle discharge, (NL)EIS, MDP) in COMSOL for a
% number of aged cell models.
%
% -- Changelog --
% 2023.09.12 | Update for gen2 toolkit | Wes H
% 2023.06.11 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath('..');
TB.addpaths;

% Constants.
cellFile = fullfile('agemodel','cellLMO_AgeSeries.mat');
TdegC = 25;
Verbose = true;
ArgHalfCycle = {0.1,100,5,TdegC,'Verbose',Verbose}; % Iavg,soc0Pct,socfPct[,...]
ArgEIS = {logspace(-3,5,50),[100 95 50 10 5],TdegC,'I',0.03,'Verbose',Verbose}; % freq,socPct[,...]
cellData = load(cellFile);

% Run simulations in COMSOL at each cell age.
clear ageSeries;
for k = length(cellData.models):-1:1
    if Verbose
        fprintf(['\nRunning RPT @ age=%.3f\n' ...
            '====================\n'],cellData.ageVect(k));
    end
    ageSeries(k) = simRPT( ...
        convertCellModel(cellData.models(k),'RLWRM'),ArgHalfCycle,ArgEIS);
end
simData.ageSeries = ageSeries;
simData.TdegC = TdegC;
simData.cellData = cellData;

% Save results to disk.
[~,modelName,~] = fileparts(cellFile);
fileName = fullfile( ...
    'simdata', ...
    'fullrpt', ...
    sprintf('%s',modelName) ...
);
save(fileName,"simData");