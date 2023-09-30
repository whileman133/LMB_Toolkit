% makeROM.m
%
% Make reduced-order model for LMB battery cell.
%
% -- Changelog --
% 2023.09.28 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath('..');
TB.addpaths;

cellFile = 'cellLMO-P2DM';
xraConfigFile = 'defaultHRA';
outdir = 'ROM_FILES';

% Generate ROM
P2DM = loadCellModel([cellFile '.xlsx']);
LLPM = convertCellModel(P2DM,'LLPM');
xraData = loadXRA([xraConfigFile '.xlsx']);
ROM = genROM(LLPM,xraData,'HRA');

% Save ROM.
if ~isfolder(outdir)
    mkdir(outdir);
end
save( ...
    fullfile(outdir,[cellFile '_' xraConfigFile '.mat']), ...
    'ROM','P2DM','LLPM' ...
);


