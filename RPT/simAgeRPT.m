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
cellFile = fullfile('agemodel','cellLMO_AgeArray.mat');
TdegC = 25;
Verbose = true;
ArgHalfCycle = {0.1,100,20,TdegC,'Verbose',Verbose}; % Iavg,soc0Pct,socfPct[,...]
ArgEIS = {logspace(-3,5,50),[100 95 50 10 5],TdegC,'I',0.03,'Verbose',Verbose}; % freq,socPct[,...]
cellData = load(cellFile);

% Run simulations in COMSOL at each cell age.
clear ageSeries; 
for i = size(cellData.ageArray,1):-1:1
    for j = size(cellData.ageArray,2):-1:1
        p = cellData.ageArray(i,j);
        if Verbose
            fprintf(['\nRunning RPT @ lam=%.3f lli=%.3f\n' ...
                '====================\n'],p.lam,p.lli);
        end
        ageSeries(i,j) = simRPT(p.model,ArgHalfCycle,ArgEIS);
    end % for
end % for
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