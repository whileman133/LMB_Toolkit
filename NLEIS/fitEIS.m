% fitEIS.m
%
% Regress linear EIS model to spectra collected in the laboratory.
%
% -- Changelog --
% 2023.06.29 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
ocpExpirementName = 'SionFresh_0C01';         % file w/ regressed OCP data
eisExpirementName = '(NL)EIS-SionCell395534'; % directory w/ raw EIS data

% Load lab impedance spectra.
spectra = loadLabNLEIS( ...
    fullfile('labdata',eisExpirementName));

% Load lab OCP.
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle');
ocpData = load( ...
    fullfile(TB.const.OCPROOT,'labdata','fitstruct',ocpExpirementName));
[TdegC,indTemp] = max(ocpData.study.testTemperatures);
ocpmodel = ocpData.study.tests(indTemp).MSMR;
ocpData = load( ...
    fullfile(TB.const.OCPROOT,'labdata','builtstruct',ocpExpirementName));
indTemp = find(ocpData.study.testTemperatures==TdegC,1,'first');
ocptest = ocpData.study.tests(indTemp);

% Perform regression.
fitLinearEIS(spectra,ocptest,ocpmodel);
