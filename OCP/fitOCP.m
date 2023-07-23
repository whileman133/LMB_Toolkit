% fitOCP.m
%
% Fit MSMR model(s) to laboratory-derived OCP estimates.
%
% -- Changelog --
% 2023.07.23 | Update for gen2 toolkit | Wesley Hileman
% 2022.08.10 | Use packaged OCP | Wesley Hileman
% 2023.07.05 | Opt. to fit one MSMR model to multiple lab | Wesley Hileman
% 2022.03.25 | Created | Wesley Hileman <whileman@uccs.edu>

clearvars -except study; clc; close all;

% Name of the study to load.
STUDYNAME = 'SionFresh_0C01';

% Directory in which to place fit models.
OUTDIR = fullfile(com.const.OCPROOT, 'labdata', 'fit');

% Directory in which to place output plots.
PLOTDIR = fullfile('..', 'plots');

% Size of the voltage bins used for computing differential capacity.
% (for more information, see: help smoothdiff.m)
VBINSIZE = 5e-3;

% Voltage range for evalulating the fit MSMR model. Should be a superset of
% the lab data collection range.
VMIN = 2; VMAX = 6;

% Number of MSMR galleries to fit.
J = 7;

% Load study data.
ocp.load(STUDYNAME,'built');

% Prepare OCP estimates for model regression.
% Augment with XPD data from the literatue.
load(fullfile(TB.const.OCPROOT,'labdata','UocpXPD.mat'));
idx = x<0.18; zXPD = flip(x(idx)); UocpXPD = flip(Uocp(idx));
estimates = ocp.Estimate.empty;
cursor = 1;
for idxTest = 1:length(builtstudy.tests)
    % Create diagonal estimate.
    test = builtstudy.tests(idxTest);
    diagEst = ocp.DiagonalEstimate(test,VBINSIZE,'maxZ',0.7,'minV',3.7);
    est = diagEst.useChg;

    % Augment with XPD data.
    Zaug = [zXPD(:); 0.18+est.Z(:)*(1-0.18)]; 
    Zaug = (Zaug-min(Zaug))/(1-min(Zaug));
    dv = est.V(1)-UocpXPD(end);
    Vaug = [UocpXPD(:)+dv; est.V(:)];
    Zunif = linspace(0,1,ceil(length(est.Z)/(1-0.18)));
    Vunif = interp1(Zaug,Vaug,Zunif,'linear','extrap');
    estAug = ocp.Estimate(test,VBINSIZE,Zunif,Vunif, ...
        sprintf('%s (Augmented with XPD Data)',est.meta));

    estimates(cursor) = estAug;
    cursor = cursor + 1;
end

% Define optimization constraints.
% - Fixed parameters: define fixed parameters using the `fix` structure.
% - Lower and upper bounds: set using either `eps` OR `lb` and `ub`
%   structures:
%       * eps (epsilon) defines symmetric lower and upper
%         bounds for a variable around the initial guess: if eps.x = a 
%         and the initial guess is x=x0, then x0*(1-a) <= x <= x0*(1+a).
%         In this way, a*100 is the percent allowed variance from the
%         initial guess.
%       * use lb and ub to specify asymmetric or arbitrary lower and upper
%         bounds, respectively.
% - Usep sets the minimum required separation of the gallery standard
%   potentials Uj0.
% - w is the relative weighting between Uocp residuals versus dtheta/dUocp
%   residuals.
opts.fix.zmin = min(zXPD);
opts.fix.zmax = 0.9999;
opts.lb.Uj0 = 3;     opts.ub.Uj0 = 5;
opts.lb.Xj = 0.02;   opts.ub.Xj = 0.8;
opts.lb.Wj = 0.01;   opts.ub.Wj = 30;
opts.Usep = 0.005;
opts.w = 0.1;

% Optimizer options. We use genetic algorithm (ga) followed by fmincon.
opts.gaPopulationSize = 200;
opts.gaIterations = 500;
opts.fminconIterations = 5000;

% Place less emphasis on XPD data in regression since it technically only
% applies at T=25degC.
opts.weighting(1).getInterval = @(U)U>4.3;
opts.weighting(1).multiplier = 0.5;
% Place more emphasis on low-potential peaks in differential capacity.
opts.weighting(2).getInterval = @(U)(3.5<U)&(U<3.7);
opts.weighting(2).multiplier = 2;

% Fit models to each OCP estimate.
models = ocp.MSMRFit.empty;
for idxEstimate = 1:length(estimates)
    estimate = estimates(idxEstimate);
    test = estimate.ocptest;

    % Take best of several optimization attempts.
    bestFit = []; 
    bestCost = Inf;
    for k = 1:10
        [fit, cost] = ocp.MSMR.fit(estimate,J,VMIN,VMAX, ...
            'fix',opts.fix,'lb',opts.lb,'ub',opts.ub,'Usep',opts.Usep,'w',opts.w, ...
            'gaPopulationSize',opts.gaPopulationSize,'gaIterations',opts.gaIterations, ...
            'fminconIterations',opts.fminconIterations,'weighting',opts.weighting);
        if cost < bestCost
            bestFit = fit;
            bestCost = cost;
        end
    end
    models(idxEstimate) = ocp.MSMRFit(test,estimate,bestFit,opts);

    ocp.MSMR.compare( ...
        3,5,2,test.temp, ...
        sprintf('Lab (%s %.0f\\circC)',test.name,test.temp),estimate, ...
        'Model Fit',bestFit);
    thesisFormat([0.2 0.1 0.1 0.1]);
    bestFit.print('Model Fit (Diagonal Estimate)');
end

fitstudy = ocp.FitStudy(builtstudy.name, models);
save(fullfile(OUTDIR,builtstudy.name),'fitstudy');