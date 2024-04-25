% compareNMC811Ds.m
clear; close all; clc;
addpath('..');
TB.addpaths;

% Load Ds vs OCP data from literature reference.
RefData = readtable('NMC811-Ds-RuessEtAl.xlsx');
UocpREF = RefData.Uocp;
DsREF_cm2secN1 = RefData.Ds;

% Load lab OCP.
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle');
ocpData = load( ...
    fullfile('..','CH6','ocp','labdata','fitstruct',ocpExpirementName));
[TdegC,indTemp] = max(ocpData.study.testTemperatures); % use highest temp.
ocpfit = ocpData.study.tests(indTemp);
warning('on','MATLAB:dispatcher:UnresolvedFunctionHandle');

figure;
semilogy(UocpREF,DsREF_cm2secN1,'o-');
thesisFormat;