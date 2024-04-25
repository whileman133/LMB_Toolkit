% makeROMFromEIS.m
%
% Make reduced-order model for LMB battery cell from parameter estimates
% derived from linear EIS.
%
% -- Changelog --
% 2023.10.05 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath('..');
TB.addpaths;

eisFile = 'EIS-16degC26degC-Ds=linear-k0=linear';
xraConfigFile = 'defaultHRA';
outdir = 'ROM_FILES';
TdegC = 25;

eisData = load(fullfile('..','eis','labfitdata',[eisFile '.mat']));
values = fastopt.evaltemp(eisData.values,eisData.modelspec,TdegC);
RLWRM = setCellParam(eisData.initialModel,values);
LLPM = convertCellModel(RLWRM,'LLPM');

% Generate ROM.
xraData = loadXRA([xraConfigFile '.xlsx']);
ROM = genROM(LLPM,xraData,'HRA');

% Save ROM.
if ~isfolder(outdir)
    mkdir(outdir);
end
save( ...
    fullfile(outdir,[eisFile '_' xraConfigFile '.mat']), ...
    'ROM','eisData','RLWRM','LLPM' ...
);