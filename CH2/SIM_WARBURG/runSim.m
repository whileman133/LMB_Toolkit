% runSim.m
%
% Run simulation study for Warburg resistance and time-constant parameters.
%
% -- Changelog --
% 2023.06.23 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath(fullfile('..','..'));
TB.addpaths('gen2');

modelName = 'cellLMO-P2DM';
p2dm = loadCellModel(modelName);
worm = convertCellModel(p2dm,'WORM');
Q = getCellParam(worm,'const.Q');

time = 0:0.1:60*10;        % time vector [sec]
iapp = zeros(size(time));  % applied current vector [A]
iapp(time>=1) = 1*Q;

modelCOMSOL = genFOM(worm);