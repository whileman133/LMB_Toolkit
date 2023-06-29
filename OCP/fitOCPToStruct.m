% fitOCPToStruct.m
%
% Utility to convert regressed OCP data to plain MATLAB structures.
%
% -- Changelog --
% 2023.05.31 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath('..');
ocp.load('SionFresh_0C01','fit');

% Convert to plain MATLAB structure and save to disk.
study = fitstudy.toStruct();
studyFile = fullfile( ...
    TB.const.OCPROOT, 'labdata', 'fitstruct', sprintf('%s.mat', fitstudy.name));
save(studyFile,'study');