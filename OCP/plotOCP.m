% plotOCP.m
%
% Plot regressed MSMR OCP curve.
%
% -- Changelog --
% 07.23.2023 | Update for gen2 toolkit | Wesley Hileman
% 10.12.2022 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath('..');
TB.addpaths();
load(fullfile(TB.const.OCPROOT,'labdata','fitstruct','FinalFit-SionFresh_0C01.mat'));
test = study.tests(study.testTemperatures==40);
electrode = MSMR(test.MSMR,'name',test.name);
electrode.plotOCP('vmin',3.49,'vmax',5);
