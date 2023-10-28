% simGITTPulse.m
%
% Simulate relaxation of cell after application of a current pulse, both
% constant-current and quarter-sine, for single step of a GITT experiment.
%
% -- Changelog --
% 2023.03.10 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
load(fullfile('simdata','GITTPulse.mat'));
Tend = ccData.Tend;
[~,indCC] = min(abs(ccData.time-Tend));
[~,indQC] = min(abs(qcData.time-Tend));

figure;
plot(ccData.time,ccData.vcell); hold on;
plot(qcData.time,qcData.vcell);
thesisFormat;

figure;
plot(ccData.output.xLocs.Thetass,ccData.output.Thetass(indCC,:)); hold on;
plot(qcData.output.xLocs.Thetass,qcData.output.Thetass(indQC,:));
thesisFormat;