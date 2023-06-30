% fitEIS.m
%
% Regress linear EIS model to spectra collected in the laboratory.
%
% -- Changelog --
% 2023.06.29 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath('..');
TB.addpaths();

% Constants.
eisExpirementName = '(NL)EIS-SionCell395534'; % directory w/ raw EIS data
ocpExpirementName = 'SionFresh_0C01';         % file w/ regressed OCP data
initialCellModelName = 'cellSionGuess-P2DM';  % model w/ initial param values

% Load lab impedance spectra.
spectra = loadLabNLEIS( ...
    fullfile('labdata',eisExpirementName));

% Load lab OCP.
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle');
ocpData = load( ...
    fullfile(TB.const.OCPROOT,'labdata','fitstruct',ocpExpirementName));
[TdegC,indTemp] = max(ocpData.study.testTemperatures);
ocpfit = ocpData.study.tests(indTemp);

% Load initial model.
initialModel = loadCellModel(initialCellModelName);

% Perform regression.
fitLinearEIS(spectra,ocpfit,initialModel);