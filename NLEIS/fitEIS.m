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
eisExpirementName = '(NL)EIS-SionCell395534_25degC'; % directory w/ raw EIS data
ocpExpirementName = 'FinalFit-SionFresh_0C01';       % file w/ regressed OCP data
initialCellModelName = 'cellSionGuess-P2DM';         % model w/ initial param values
solidDiffusionModel = 'msmr';  % selects solid diffusion model for porous electrode
kineticsModel = 'linear';      % selects kinetics model for porous electrode

% Load lab impedance spectra.
spectra = loadLabNLEIS( ...
    fullfile('labdata',eisExpirementName));

% Load lab OCP.
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle');
ocpData = load( ...
    fullfile(TB.const.OCPROOT,'labdata','fitstruct',ocpExpirementName));
[TdegC,indTemp] = max(ocpData.study.testTemperatures); % use highest temp.
ocpfit = ocpData.study.tests(indTemp);
warning('on','MATLAB:dispatcher:UnresolvedFunctionHandle');

% Load initial model.
initialModel = loadCellModel(initialCellModelName);

% Perform regression.
fitData = fitLinearEIS(spectra,ocpfit,initialModel, ...
    'SolidDiffusionModel',solidDiffusionModel, ...
    'KineticsModel',kineticsModel,...
    'WeightFcn',@getWeight);

% Save results to disk.
fileName = fullfile( ...
    'labfitdata', ...
    sprintf('EIS-Cell%s-%.0fdegC-Ds=%s-k0=%s', ...
            spectra.cellName,spectra.TdegC,solidDiffusionModel,kineticsModel) ...
);
save(fileName,'-struct',"fitData");

% Residual weighting function. Returns relative weight corresponding to the
% specified frequency and SOC setpoint.
function w = getWeight(freq,socPct)
    w = 1;
%     if socPct<=25
%         % Linear derate from 25% to 0% SOC.
%         w = socPct/25;
%     else
%         w = 1;
%     end
end

