% saveParamValues.m
%
% Collect and save parameter values determined by TF model regression.
%
% -- Changelog --
% 2024.03.08 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath(fullfile('..','..'));
TB.addpaths;

% Name of file to process.
file = fullfile('labfitdata','EIS-16degC26degC-Ds=linear-k0=linear.mat');

% Reference temperature.
TrefdegC = 25;

soc = linspace(0,1,100);
socPct = soc*100;
fitData = load(file);
TdegC = fitData.TdegC;
values = fastopt.unpack( ...
    fastopt.pack(fitData.values,fitData.modelspec),fitData.modelspec, ...
    'flat',true,'sparse',true);
[values, EactStruct] = fastopt.evaltemp(values,fitData.modelspec,TrefdegC,true);
lb = fastopt.unpack( ...
    fastopt.pack(fitData.lb,fitData.modelspec),fitData.modelspec, ...
    'flat',true,'sparse',true);
ub = fastopt.unpack( ...
    fastopt.pack(fitData.ub,fitData.modelspec),fitData.modelspec, ...
    'flat',true,'sparse',true);
cellName = fitData.arg.labSpectra.cellName;
cellName = strsplit(cellName,'_');
cellName = cellName{1};
labels = sprintf('T=%.0f \\circC (#%s)',TdegC,cellName);

VarName = fieldnames(values);
Value = cellfun(@(x)values.(x)(1),VarName,'UniformOutput',false);
LB = cellfun(@(x)lb.(x)(1),VarName,'UniformOutput',false);
UB = cellfun(@(x)ub.(x)(1),VarName,'UniformOutput',false);
EactKJ = cellfun(@(x)EactStruct.(x)/1000,VarName,'UniformOutput',false);
tab = table(VarName,Value,LB,UB,EactKJ);
writetable( ...
    tab,fullfile('labfitdata-xlsx','LinearEIS-AllCells.xlsx'), ...
    "Sheet",'values');