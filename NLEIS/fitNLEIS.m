% fitEIS.m
%
% Regress nonlinear EIS model to spectra collected in the laboratory.
%
% -- Changelog --
% 2023.09.23 | Fix SOC computation error | Wes H.
% 2023.08.31 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath('..');
TB.addpaths();

% Constants.
linearFitName = '202310_EIS-16degC26degC-Ds=linear-k0=linear';
eisExpirementNames = {
    '(NL)EIS-SionCell395534_15degC'  % directories w/ raw EIS data
    '(NL)EIS-SionCell395534_25degC'
};
useParallel = true;             % enable/disable parallel processing

% Load lab impedance spectra.
clear spectra;
for m = length(eisExpirementNames):-1:1
    spectra(m) = loadLabNLEIS(fullfile('labdata',eisExpirementNames{m}));
end

% Load linear fit data.
labLinearFit = load(fullfile('labfitdata',[linearFitName '.mat']));

% Perform regression.
fitData = fitNonlinearEIS(spectra,labLinearFit, ...
    'WeightFcn',@getWeight,...
    'UseParallel',useParallel);

% Save results to disk.
tempstr = sprintf('%.0fdegC',[spectra.TdegC]);
fileName = fullfile( ...
    'labfitdata', ...
    sprintf('NLEIS-%s',tempstr) ...
);
save(fileName,'-struct','fitData');

% Residual weighting function. Returns relative weight corresponding to the
% specified frequency, SOC setpoint, and temperature.
function w = getWeight(freq,socPct,TdegC)
    w = 1;
    if freq<=10e-3
        w = 2; % compensate for lower point density below 10mHz
    end
    if freq>10
        w = 0; % examination of phase indicates mostly noise above 10Hz;
               % do not include these data-points in model regression
    end
end