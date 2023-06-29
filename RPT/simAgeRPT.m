% simAgeRPT.m
%
% Simulate full RPT (half-cycle discharge, (NL)EIS, MDP) in COMSOL for a
% number of aged cell models.
%
% -- Changelog --
% 2023.06.11 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
if ~exist('TOOLBOX_LMB','dir')
    % bootstrap the toolbox
    addpath('..');
    TB.addPaths();
end

% Constants.
cellFile = 'cellLMO_AgeSeries.mat';
TdegC = 25;
Verbose = true;
ArgHalfCycle = {0.1,100,5,TdegC,'Verbose',Verbose}; % Iavg,soc0Pct,socfPct[,...]
ArgEIS = {logspace(-3,5,50),[100 95 50 10 5],TdegC,'I',0.03,'Verbose',Verbose}; % freq,socPct[,...]
ArgPulse = {0.01,30,50,TdegC,'Verbose',Verbose}; % I,Tp,socPct[,...]
cellData = load(cellFile);

% Run simulations in COMSOL at each cell age.
clear ageSeries;
for k = length(cellData.lumped):-1:1
    if Verbose
        fprintf(['\nRunning RPT @ age=%.3f\n' ...
            '====================\n'],cellData.ageVect(k));
    end
    ageSeries(k) = simRPT( ...
        cellData.lumped(k),ArgHalfCycle,ArgEIS,ArgPulse);
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
if ~isempty(suffix)
    fileName = [fileName '-' suffix];
end
save(fileName,"simData");