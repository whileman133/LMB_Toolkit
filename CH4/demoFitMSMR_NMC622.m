% testFitMSMR_NMC622.m
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
truth = MSMR.NMC622();  % true MSMR model of electrode
J = truth.J;   % Number of galleries to fit (assume foud by trial-and-error)
w = 0.1;       % Weight of differential capacity in cost function [-]
Usep = 0.005;  % Minimum required separation of gallery potentials U0 [V]
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
gaPopulationSize = 200;
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
lb.U0 = 3;         ub.U0 = 5;
lb.X  = 0.1;       ub.X  = 0.4;
lb.omega = 0.1;    ub.omega = 10;
lb.thetamin = 0;   ub.thetamin = 0.2;
lb.thetamax = 0.9; ub.thetamax = 1.0;


% Regress MSMR model ------------------------------------------------------

data = fitMSMR(ocp,J, ...
    'lb',lb,'ub',ub,'w',w,'Usep',Usep, ...
    'gaPopulationSize',gaPopulationSize, ...
    'gaIterations',gaIterations, ...
    'fminconIterations',fminconIterations, ...
    'verbose',verbose);

% Print results.
printMSMR('Fit Model',data.est);
printMSMR('Truth',truth.toStruct());
printMSMR('Percent Error',getPercentError(data.est,truth.toStruct()));

% Save results to disk.
save('demoFitMSMR_NMC622.mat','data','truth');


function printMSMR(name,params)
%PRINTMSMR Print MSMR parameters to the console.

fprintf("%s\n",upper(name));
fprintf("  U0      : ");
fprintf("%-7.3f",params.U0);
fprintf("\n");
fprintf("  X       : ");
fprintf("%-7.3f",params.X);
fprintf("\n");
fprintf("  omega   : ");
fprintf("%-7.3f",params.omega);
fprintf("\n");
fprintf("  thetamin: ");
fprintf("%-7.3f",params.thetamin);
fprintf("\n");
fprintf("  thetamax: ");
fprintf("%-7.3f",params.thetamax);
fprintf("\n");

end

function [err, rmse] = getPercentError(estimate,truth)
%GETRMSE Calculate RMSE between estimates and true parameter values.

paramnames = fieldnames(truth);
toterr = 0;
totvar = 0;
for k = 1:length(paramnames)
    pname = paramnames{k};
    tru = truth.(pname);
    est = estimate.(pname);
    err.(pname) = 100*(tru-est)./tru;
    toterr = toterr + sum(err.(pname).^2);
    totvar = totvar + length(err.(pname));
end
rmse = sqrt(toterr/totvar);

end