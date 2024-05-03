% compareNMC811Ds.m
clear; close all; clc;
addpath('..');
TB.addpaths;
ocpExpirementName = 'Fit-SionFresh_0C01';  % file w/ regressed OCP data
gprEstimateName = 'EIS-16degC26degC-Ds=linear-k0=linear_DsEstimate';
TdegC = 25;

% Load Ds vs OCP data from literature reference.
dataREF = readtable('NMC811-Ds-RuessEtAl.xlsx');
UocpREF = dataREF.Uocp;
DsREF = dataREF.Ds/100/100; % convert units from cm^2/s to m^2/s

% Load lab OCP.
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle');
ocpData = load( ...
    fullfile('..','CH6','ocp','labdata','fitstruct',ocpExpirementName));
[~,indTemp] = max(ocpData.study.testTemperatures); % use highest temp.
ocpfit = ocpData.study.tests(indTemp);
paramsOCP = ocpfit.MSMR;
warning('on','MATLAB:dispatcher:UnresolvedFunctionHandle');

% Load lab Ds estimate.
dataEIS = load(fullfile('..','CH6','eis','gpr',gprEstimateName));
thetaEIS = dataEIS.est.theta;
DsEIS = 10.^dataEIS.est.log10Ds;
sigEIS = 10.^dataEIS.est.sigma;

% Load GITT data.
gittData = load('GITT-1.mat');
thetaGITT = gittData.est.theta;
DsGITT = 10.^gittData.est.log10Ds;

% Move DsREF over stiochometry.
ocpREF = MSMR(paramsOCP).ocp('voltage',UocpREF,'TdegC',TdegC);
thetaREF = ocpREF.theta;

% Estimate particle radius as scaling factor.
% Note: Ds = Ds(bar)*Rs^2.
DsEIS0 = interp1(thetaEIS,DsEIS,thetaREF); % move EIS Ds over ref. theta.
% !! Perform the LS regression over logspace
y = log10(DsREF./DsEIS0);
M = ones(size(DsREF));
Rs = 10.^((M\y)/2);  % M\y is equivalent to mean(y) in this case

% Remove scaling factor from GITT estimate.
% Note: Ds = Ds(bar)*k; 
DsGITT0 = interp1(thetaGITT,DsGITT,thetaREF); % move GITT Ds over ref. theta.
y = log10(DsREF./DsGITT0);
M = ones(size(DsREF));
k = 10.^(M\y);  % M\y is equivalent to mean(y) in this case

% Plot result.
figure;
semilogy(thetaEIS,DsEIS*Rs^2); hold on;
semilogy(thetaGITT,DsGITT*k,':');
semilogy(thetaREF,DsREF,'o-');
set(gca,'xdir','reverse');
xlabel('Li Composition, \theta_s [-]');
ylabel('Solid Diffusivity, D_s [m^2s^{-1}]');
title('Solid Diffusivity of NMC811 vs. Composition');
legend( ...
    sprintf('EIS (R_s \\approx %.2f\\mum)',1e6*Rs), ...
    sprintf('GITT'), ...
    sprintf('Ruess \\itet al.\\rm (GITT) '), ...
    'Location','best');
thesisFormat;