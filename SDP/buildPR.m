% makePR.m
%
% Pre-process raw pulse-resistance (PR) data collected by a Gamry 
% potentiostat into a format usable by other scripts.
%
% Type `help pr.build` for more information.
%
% -- Changelog --
% 08.16.2022 | Improved packaging, naming, documentation | Wesley Hileman
% 04.20.2022 | Created
% Wesley Hileman <whileman@uccs.edu>
% University of Colorado Colorado Springs

clear; clc; close all;

% Name of the study to build.
STUDYNAME = 'SionFresh';

% Pulse-current setpoints (not exact/measured, but what was programmed) [A]
% Note: these should be positive (symmetric negative values assumed) and
% should appear in the order traversed by the Gamry script.
CURRENT = 0.01:0.03:0.28;

% Timing specs.
T0 = 0;       % Time at onset of the pulse [s].
TF = 0.1;     % Time at offset of the pulse [s].
FS = 300000;  % Nominal sampling rate [Hz].

% Process files on disk.
study = pr.build(STUDYNAME, CURRENT, T0, TF, FS); %...
    %'only', {'Cell395535'}, 'merge', true);