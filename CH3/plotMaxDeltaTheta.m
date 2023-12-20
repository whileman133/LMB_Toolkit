% plotMaxDeltaTheta.m
%
% Plot the maximum composition and SOC deviation allowed such that 
% error in the linear approximation of the OCP is limited to some DeltaU.
%
% -- Changelog --
% 2023.12.19 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath('..'); TB.addpaths;

% Constants.
deltaU = 1e-3;       % [V]
deltaThetaB = 0.05;  % [-]
electrode = MSMR.NMC622;
theta0 = electrode.zmax;
theta100 = electrode.zmin;


% Computations ------------------------------------------------------------
ocpData = electrode.ocp();
theta = ocpData.theta;
soc = (theta-theta0)./(theta100-theta0);
Uocp = ocpData.Uocp;
d2Uocp = ocpData.d2Uocp;
deltaThetaMax = min(deltaThetaB,sqrt(2*deltaU./abs(d2Uocp)));
deltaSocMax = deltaThetaMax/abs(theta100-theta0);


% Plotting ----------------------------------------------------------------
figure;
plot(theta,deltaThetaMax);
xlabel('Fractional Composition, \theta_s');
ylabel('Max. Compisition Delta, \Delta\theta_{max}');
title('Maximum Composition Deviation (NMC622)');
thesisFormat;
print(fullfile('plots','maxDeltaTheta'),'-depsc');
print(fullfile('plots','maxDeltaTheta'),'-dpng');

% figure;
% plot(theta,Uocp);
% thesisFormat;