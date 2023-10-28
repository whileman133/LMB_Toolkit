% fitEIS.m
%
% Regress linear EIS model to spectra collected in the laboratory.
%
% NOTE: With useParallel=false, this regression saves multiple solutions
% from the optimization. This is needed for estimating diffusivity Ds and
% charge-transfer resistance in the runGPRLinearEIS.m script.
%
% -- Changelog --
% 2023.09.23 | Fix SOC computation error | Wes H.
% 2023.06.29 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath('..');
TB.addpaths();

% Constants.
eisExpirementNames = {
    '(NL)EIS-SionCell395534_15degC'  % directories w/ raw EIS data
    '(NL)EIS-SionCell395534_25degC'
};
ocpExpirementName = 'FinalFit-SionFresh_0C01';  % file w/ regressed OCP data
initialCellModelName = 'cellSionGuess-P2DM';    % model w/ initial param values
solidDiffusionModel = 'linear';  % selects solid diffusion model for porous electrode
kineticsModel = 'linear';        % selects kinetics model for porous electrode
useParallel = false;              % enable/disable parallel processing
prefix = '202310_';

% Load lab impedance spectra.
clear spectra;
for m = length(eisExpirementNames):-1:1
    spectra(m) = loadLabNLEIS(fullfile('labdata',eisExpirementNames{m}));
end

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
    'WeightFcn',@getWeight,...
    'UseParallel',useParallel);

% Save results to disk.
tempstr = sprintf('%.0fdegC',[spectra.TdegC]);
fileName = sprintf( ...
    'EIS-%s-Ds=%s-k0=%s',tempstr,solidDiffusionModel,kineticsModel);
if ~isempty(prefix)
    fileName = [prefix fileName];
end
save(fullfile('labfitdata',fileName),'-struct','fitData');

% Residual weighting function. Returns relative weight corresponding to the
% specified frequency, SOC setpoint, and temperature.
function w = getWeight(freq,socPct,TdegC)
    w = 1;
    if freq<=10e-3
        w = 2; % compensate for lower point density below 10mHz
    end
end