% plotPRSurfaceTheoretical.m
%
% Compute and plot a theoretical pulse-resistance surface given parameter 
% values of a hypothetical model. Useful for getting a feel of the pulse
% resistance surface.
%
% -- Changelog --
% 07.23.2023 | Updated for gen2 toolkit | Wesley Hileman
% 08.16.2022 | Updated for new schema, use real OCP params | Wesley Hileman
% 03.18.2022 | Created |
% Wesley Hileman <whileman@uccs.edu>
% University of Colorado Colorado Springs

clear; clc; close all;
addpath('..');
TB.addpaths;
addpath(TB.const.OCPROOT);

% Define SOC/current setpoints.
soc = linspace(0.01, 0.95, 20);
iapp = linspace(-0.3, 0.3, 20);
TdegC = 25;
T = TdegC+273.15;

% Construct pulse model.
model = loadCellModel('cellSionGuess-P2DM.xlsx');
model = convertCellModel(model,'LLPM'); % legacy lumped parameter model
msmr = MSMR(model.function.pos);

% Compute pulse resistance and associated quantities.
[R0, deltaV, phi_se0, phi_se2, phi_se3] = ...
    pr.getPulseResistance(model,soc,iapp,TdegC);

% Plot results.
figure; surf(soc*100,iapp,R0); 
xlabel('State of Charge [%]'); 
ylabel('Pulse Magnitude [A]');
zlabel('Pulse Resistance, R_0 [\Omega]');
thesisFormat3d;
