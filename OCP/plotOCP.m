% plotOCP.m

clear; close all; clc;
addpath('..');
TB.addpaths();
load(fullfile(TB.const.OCPROOT,'labdata','fitstruct','FinalFit-SionFresh_0C01.mat'));
test = study.tests(study.testTemperatures==40);
electrode = MSMR(test.MSMR,'name',test.name);
electrode.plotOCP('vmin',3.49,'vmax',5);
