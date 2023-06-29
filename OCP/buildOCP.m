% buildOCP.m
%
% Pre-process raw OCP data collected by a Gamry potentiostat into a 
% format usable by other scripts.
% 
% Type `help ocp.build` for more information.
%
% -- Changelog --
% 08.10.2022 | Move OCP data directory, various cleanup | Wesley Hileman
% 08.04.2022 | Collect into package, rename to buildOCP | Wesley Hileman
% 06.29.2022 | Parse a single study at a time | Wesley Hileman
% 02.10.2022 | Created |
% Wesley Hileman <whileman@uccs.edu>
% University of Colorado Colorado Springs

clear; clc; close all;
addpath('..');

% Name of the study to build.
STUDYNAME = 'SionFresh_0C01';

% Reference temprature for calibration scripts [C]. Ensure there is at
% least one test performed at the reference temperature.
TREF = 25;

% Cell voltage limits [V]. Not used in processing, for later reference
% purposes only.
VMIN = 3.2;
VMAX = 4.3;

% Process files on disk.
built = ocp.build(STUDYNAME, TREF, VMIN, VMAX, 'verbose', true);

% Convert to plain MATLAB structure and save to disk.
study = built.toStruct();
studyFile = fullfile( ...
    TB.const.OCPROOT, 'labdata', 'builtstruct', sprintf('%s.mat', built.name));
save(studyFile,'study');