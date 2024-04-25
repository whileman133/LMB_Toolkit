% fitOCP.m
%
% Fit MSMR model(s) to laboratory-derived OCP estimates.
% This implementation uses ga.m (genetic algorithm) followed by fmincon.m
% to regress the MSMR model to laboratory OCP collected at several
% temperatures (T=0degC, T=15degC, T=25degC, T=40degC, T=55degC).
% **We consider the T=40degC dataset most indicative of the true OCP.**
%
% XPD data from the literature (Friedrich et al.) is also used to extend 
% the range of OCP measurements over lower lithiation than can be permitted 
% by our laboratory tests. This allows for better fitting of the MSMR 
% galleries.
%
% -- Changelog --
% 2023.07.23 | Update for gen2 toolkit | Wesley Hileman
% 2022.08.10 | Use packaged OCP | Wesley Hileman
% 2023.07.05 | Opt. to fit one MSMR model to multiple lab | Wesley Hileman
% 2022.03.25 | Created | Wesley Hileman <whileman@uccs.edu>

clear; clc; close all;
addpath(fullfile('..','..'));
TB.addpaths;

% Name of the study to load.
STUDYNAME = 'SionFresh_0C01';  % C/100 quasi-equilibrium cycling

% Directory in which to place fit models.
OUTDIR = fullfile(TB.const.OCPROOT, 'labdata', 'fit');
OUTNAME = ['Fit-' STUDYNAME];

% Directory in which to place output plots.
PLOTDIR = fullfile('.', 'plots', ['FIT_' STUDYNAME]);
if ~isfolder(PLOTDIR), mkdir(PLOTDIR); end

% Size of the voltage bins used for computing differential capacity.
% (for more information, see: help smoothdiff.m)
VBINSIZE = 8e-3;

% Voltage range for evalulating the fit MSMR model. Should be a superset of
% the lab data collection range.
VMIN = 2; VMAX = 6;

% Number of MSMR galleries to fit.
J = 7;

% Number of times to run the ga/fmincon optimization at each temperature.
% The run with the best cost will be used as the final result.
NUMBER_OF_ATTEMPTS = 1;

% Temperatures to process.
TDEGC = [40];

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

    if ~any(test.temp == TDEGC)
        continue;
    end

    % Plot differential Capacity.
    figure;
    plot(diagEst.unifZ,diagEst.disdzdvUnifZ,'r'); hold on;
    plot(diagEst.unifZ,diagEst.chgdzdvUnifZ,'b');
    xlabel('SOC');
    ylabel('AbsDiffCap');
    title('DiffCap-v-SOC');
    legend('Discharge','Charge','Location','southeast');
    thesisFormat;
    print(fullfile(PLOTDIR,sprintf('DiffCapDIS-%.0fdegC',test.temp)),'-depsc');
    print(fullfile(PLOTDIR,sprintf('DiffCapDIS-%.0fdegC',test.temp)),'-dpng');
    figure;
    plot(diagEst.disdzdvUnifV,diagEst.unifV,'r'); hold on;
    plot(diagEst.chgdzdvUnifV,diagEst.unifV,'b');
    xlabel('AbsDiffCap');
    ylabel('Voltage');
    title('DiffCap-v-SOC');
    legend('Discharge','Charge','Location','southwest');
    thesisFormat;
    print(fullfile(PLOTDIR,sprintf('DiffCapCHG-%.0fdegC',test.temp)),'-depsc');
    print(fullfile(PLOTDIR,sprintf('DiffCapCHG-%.0fdegC',test.temp)),'-dpng');

    % Plot OCP estimate.
    figure;
    ind = 3.2<=est.V&est.V<=4.3;
    plot(test.disZ,test.disV,'r:'); hold on;
    plot(test.chgZ,test.chgV,'b:');
    plot(est.Z(ind),est.V(ind),'k');
    ylim([3.2 4.3]);
    xlabel('SOC');
    ylabel('Voltage');
    title('OCP');
    l = legend('Vdis','Vchg','Uocp','Location','southeast');
    l.Position(1) = l.Position(1) - 0.06;
    l.FontSize = 20;
    thesisFormat;
    addInset([0 0.1],[0.07 3.35]);
    addInset([0.9 1.0],[0.6 3.85],'YSpan',[3.4 3.55]);
    print(fullfile(PLOTDIR,sprintf('UocpHAT-%.0fdegC',test.temp)),'-depsc');
    print(fullfile(PLOTDIR,sprintf('UocpHAT-%.0fdegC',test.temp)),'-dpng');

    % Augment with XPD data.
    Zaug1 = [zXPD(:); 0.18+est.Z(:)*(1-0.18)]; 
    Zaug = (Zaug1-min(Zaug1))/(1-min(Zaug1));
    dv = est.V(1)-UocpXPD(end);
    Vaug = [UocpXPD(:)+dv; est.V(:)];
    Zunif = linspace(0,1,ceil(length(est.Z)/(1-0.18)));
    Vunif = interp1(Zaug,Vaug,Zunif,'linear','extrap');
    estAug = ocp.Estimate(test,VBINSIZE,Zunif,Vunif, ...
        sprintf('%s (Augmented with XPD Data)',est.meta));

    % Plot augmented OCP estimate.
    figure;
    plot(Zaug1,Vaug);
    xline(0.18,'k');
    yline(4.3,'k');
    ylim([3.2 4.6]);
    xlabel('theta');
    ylabel('Voltage');
    title('OCP');
    thesisFormat;
    print(fullfile(PLOTDIR,sprintf('UocpHATAUG-%.0fdegC',test.temp)),'-depsc');
    print(fullfile(PLOTDIR,sprintf('UocpHATAUG-%.0fdegC',test.temp)),'-dpng');

     % Plot differential capacity of augmented OCP estimate.
    figure;
    plot(estAug.refZ,estAug.dzdvRefZ);
    xlabel('SOC');
    ylabel('AbsDiffCap');
    title('DiffCap-v-SOC');
    thesisFormat;
    print(fullfile(PLOTDIR,sprintf('DiffCapAUG-Z-%.0fdegC',test.temp)),'-depsc');
    print(fullfile(PLOTDIR,sprintf('DiffCapAUG-Z-%.0fdegC',test.temp)),'-dpng');
    figure;
    plot(estAug.dzdvRefV,estAug.refV);
    ylim([test.vmin max(estAug.V)]);
    xlabel('AbsDiffCap');
    ylabel('Voltage');
    title('DiffCap-v-SOC');
    thesisFormat;
    print(fullfile(PLOTDIR,sprintf('DiffCapAUG-U-%.0fdegC',test.temp)),'-depsc');
    print(fullfile(PLOTDIR,sprintf('DiffCapAUG-U-%.0fdegC',test.temp)),'-dpng');

    estimates(cursor) = estAug;
    cursor = cursor + 1;
end

return;

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
    for k = 1:NUMBER_OF_ATTEMPTS
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

% Save regressed models to disk.
fitstudy = ocp.FitStudy(builtstudy.name, models);
save(fullfile(OUTDIR,OUTNAME),'fitstudy');

% Convert to plain MATLAB structure and save to disk.
study = fitstudy.toStruct();
studyFile = fullfile( ...
    TB.const.OCPROOT,'labdata','fitstruct',sprintf('%s.mat',OUTNAME));
save(studyFile,'study');