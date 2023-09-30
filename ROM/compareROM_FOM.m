% compareROM_FOM.m
%
% Compare predictions of ROM to full-order COMSOL model.
%
% -- Changelog --
% 2023.09.29 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath('..');
TB.addpaths;

romFile = 'cellLMO-P2DM_defaultHRA';
romData = load(fullfile('ROM_FILES',[romFile '.mat']));
ROM = romData.ROM;
LLPM = romData.LLPM;
Q = getCellParams(LLPM,'const.Q');

% Constants.
socPct0 = 100;  % starting SOC [%]
Tstart = 10;    % time at application of current [s]
Tend = 3600;    % time at end of simulation [s]
I = Q/2;        % discharge current magnitude [A]
TdegC = 25;     % cell temperature [degC]

% Build simulation waveforms.
time = 0:1:Tend;
iapp = zeros(size(time));
iapp(time>=Tstart) = I;
Tvect = TdegC*ones(size(time));
simData.SOC0 = socPct0;
simData.Iapp = iapp;
simData.T = Tvect;
simData.time = time;
simData.TSHIFT = 0;

% Simulate ROM.
ROMout = simROM(ROM,simData,'outBlend');

% Simulate FOM.
genData = genFOM(LLPM);
[FOM,FOMout] = simFOM(genData,simData);

% Save results.
save(fullfile('SIM_FILES',[romFile '.mat']),'FOMout','ROMout','simData','romData');
