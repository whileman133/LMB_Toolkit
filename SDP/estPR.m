% estPR.m
%
% Estimate resistance from pulse-response measurements.
% The output consists of resistance surface esstimates at each temperature.
%
% -- Changelog --
% 07.23.2023 | Update for gen2 toolkit | Wesley Hileman
% 08.16.2022 | Packaging and naming improvements | Wesley Hileman
% 06.27.2022 | Upated | Wesley Hileman
% 03.24.2022 | Created |
% Wesley Hileman <whileman@uccs.edu>
% University of Colorado Colorado Springs

% Takes a long time to load pulse data, don't clear it.
clearvars -except prstudy; clc; close all;
addpath('..');
TB.addpaths;
addpath(TB.const.OCPROOT);

PRSTUDYNAME = 'SionFresh';        % PR study containing raw pulse data.
OCPSTUDYNAME = 'SionFresh_0C01';  % OCP study containing capacity estimates.

% Directory in which to place processed estimate data.
OUTDIR = fullfile(TB.const.SDPROOT,'labdata','surfaces');

% Directory in which to save plots.
PLOTDIR = fullfile(TB.const.SDPROOT,'labdata','plots','surfaces');

% Load studies.
pr.load(PRSTUDYNAME,'built','varname','prstudy');
ocp.load(OCPSTUDYNAME,'built','varname','ocpstudy');

% ECM model fit settings.
tmax = 100e-6;    % Pulse duration over which to regress model.
eps.taux = 20;    % Allowed variance on time constants from initial guess (ratiometric - multiplier/divider).
taumin = 1e-6;    % Minimum allowed time constant (needs to be a bit >0 to prevent rank deficiency warnings).
nRC = 2;          % Number of RC time constants to fit.

% Residual weighting. Define time interval(s) over which residuals are
% weighted higher than normal. (Provide a structure array for more than
% one weighting interval).
weighting(1).t1 = 0;
weighting(1).t2 = 100e-6;
weighting(1).w = 100;
weighting(2).t1 = 100e-6;
weighting(2).t2 = 2e-3;
weighting(2).w = 10;

% Upper optimization bounds.
ub.taux = ones(nRC,1);  % Defaults
% Limit lowest time constant to several us maximum, otherwise we can "miss" 
% fast dynamics.
ub.taux(1) = 50e-6; 

% Curve fit analysis ------------------------------------------------------
% T = 0; soc = 40; ipulse = -0.1;
% test = prstudy.getTest(T);
% pulses = test.getPulses(soc,ipulse);
% ncol = 5; nrow = ceil(length(pulses)/ncol);
% disptmax = tmax;
% saveplots = false;
% 
% pulse = pulses(2);
% guess = pr.ECM_xRC.guess(pulse,tmax,nRC);
% [est, J, tp, ip, vp] = pr.ECM_xRC.fit( ...
%     pulse,guess,tmax,eps,taumin,'weighting',weighting,'ub',ub);
% vg = guess.getV(ip,tp);
% vm = est.getV(ip,tp);
% maxV = max(max(vp),max(max(vg),max(vm)));
% figure;
% plot(tp*1e3,vp*1e3,'-d','MarkerSize',3,'MarkerFaceColor','auto'); hold on;
% %plot(tp*1e6,vg*1e3,':d','MarkerSize',3,'MarkerFaceColor','auto');
% plot(tp*1e3,vm*1e3,':d','MarkerSize',3,'MarkerFaceColor','auto');
% xlim([0 disptmax*1e3]); %ylim([0 400]);
% xlabel('Time, t [ms]');
% ylabel('\Deltav [mV]');
% title(sprintf('Pulse Response (i_{app}=%.0fmA, SOC=%d%%, T=%.0f\\circC)', ...
%     ipulse*1000,soc,test.temp));
% %legend('Lab','Guess','Model Fit','Location','southeast');
% legend('Lab','Model Fit','Location','southeast');
% if saveplots
%     filename = fullfile(PLOTDIR, ...
%         sprintf('PRfit-%s-%s-T%-.0f-I%.0fmA-soc%d-tmax%.0fus', ...
%             prstudy.name,test.name,test.temp,ipulse*1000,soc,disptmax*1e6) ...
%     );
%     todisk(gcf,filename);
% end

% Pulse Resistance Computation --------------------------------------------
tests = prstudy.tests;
currents = tests(1).currents(:); socs = tests(1).socs(:);
surfaces = pr.Surface.empty(0,length(tests));
for idxTest = 1:length(tests)
    prtest = tests(idxTest);
    try
        ocptest = ocpstudy.getTest(prtest.temp);
    catch
        % Default to using T=25degC OCP data.
        ocptest = ocpstudy.getTest(25);
    end
    ECMs = cell(length(currents),length(socs));

    for idxSOC = 1:length(socs)
        soc = socs(idxSOC);
        fprintf('%3d%% SOC (%3d of %3d)\n', soc, idxSOC, length(socs));
    
        for idxCurrent = 1:length(currents)
            I = currents(idxCurrent);
            pulses = prtest.getPulses(soc,I);
            est = pr.ECM_xRC.empty(0,length(pulses));
            fprintf('\t%10fA (%3d of %3d) ', I, idxCurrent, length(currents));
    
            for idxRepeat = 1:length(pulses)
                pulse = pulses(idxRepeat);
                guess = pr.ECM_xRC.guess(pulse,tmax,nRC);
                [est(idxRepeat), J, tp, ip, vp] = pr.ECM_xRC.fit( ...
                    pulse,guess,tmax,eps,taumin,'weighting',weighting,'ub',ub);
                fprintf('.');
            end
            fprintf('\n');

            ECMs{idxCurrent,idxSOC} = est;
        end
    end

    opts.tmax = tmax; 
    opts.nRC = nRC; 
    opts.taumin = taumin; 
    opts.eps = eps;
    % Need ocptest for capacity estimates to determine absolute SOC
    est = pr.Surface(prtest,ocptest,ECMs,opts);
    surfaces(idxTest) = est;
end
surfaces = pr.Surfaces(prstudy.name,surfaces);
studyFile = fullfile( ...
    OUTDIR, ...
    sprintf('%s-%dtau-%dus.%s', PRSTUDYNAME, nRC, round(tmax*1e6), 'mat'));
save(studyFile,'surfaces');
