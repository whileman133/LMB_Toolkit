% testFitMSMR_C6.m
%
% Verify the MSMR model regression using synthetic data.
%
% -- Changelog -- 
% 2024.01.21 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath('..');
TB.addpaths;
rng(0);  % make results repeatable


% Constants ---------------------------------------------------------------
TdegC = 25;             % temperature [degC]
truth = MSMR.C6();  % true MSMR model of electrode
J = truth.J;   % Number of galleries to fit (assume foud by trial-and-error)
w = 0.1;       % Weight of differential capacity in cost function [-]
Usep = 0.005;  % Minimum required separation of gallery potentials U0 [V]
tries = 1000;  % Number of independent runs of the optimization
verbose = true;


% Generate synthetic OCP --------------------------------------------------
% Evalulate OCP-vs-relative-composition curve.
ocpData = truth.ocp('npoints',1000);
ocp.TdegC = TdegC;
ocp.Z = (ocpData.theta-truth.thetamin)/(truth.thetamax-truth.thetamin);
ocp.U = ocpData.Uocp;
ocp.dZ = (1./ocpData.dUocp)/(truth.thetamax-truth.thetamin); % !!! important to scale


% Optimization control ----------------------------------------------------
% Genetic algorithm and fmincon parameters
gaPopulationSize = 300;
gaIterations = 200;
fminconIterations = 5000;

% Optimization bounds:
% * tighter bounds increase the liklihood of regression success
% * rough limits on U0, X, and omega can usually be gleamed from
%   plots of the differential capacity vs Uocp and SOC; you can use
%   the inspectMSMR.m function to help
% * bounds can be imposed on the individual galleries by making
%   U0, X, and omega column vectors; the optimization is garenteed to
%   preserve the order of the gallery potentials U0
lb.U0 = 0.05;      ub.U0 = 0.4;
lb.X  = 0.05;      ub.X  = 0.5;
lb.omega = 0.05;   ub.omega = 10;
lb.thetamin = 0;   ub.thetamin = 0.1;
lb.thetamax = 0.9; ub.thetamax = 1.0;


% Regress MSMR model ------------------------------------------------------

data = fitMSMR(ocp,J, ...
    'lb',lb,'ub',ub,'w',w,'Usep',Usep, ...
    'gaPopulationSize',gaPopulationSize, ...
    'gaIterations',gaIterations, ...
    'fminconIterations',fminconIterations, ...
    'verbose',verbose);

% Save results to disk.
save('demoFitMSMR_C6.mat','data','truth');
