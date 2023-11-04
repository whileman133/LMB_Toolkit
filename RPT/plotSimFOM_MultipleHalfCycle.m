% plotSimFOM_MultipleHalfCycle.m
%
% Plot simulated half-cycle discharge, compare to Baker and Verbrugge
% perturbation approximation.
%
% -- Changelog --
% 2023.06.01 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath('..');
TB.addpaths();

% Load simulation data.
simName = 'cellLMO-P2DM-100pct-0pct';
load(fullfile('simdata','multhalfcyc',[simName '.mat']));

% Fetch time, iapp(t), vcell(t) vectors.
halfcycData = [simData.iappSeries];
halfcycArg = [halfcycData.param];
soc0Pct = [halfcycArg.soc0Pct].';
Iapp = [halfcycData.iapp];
Vcell = [halfcycData.vcell];
time = zeros(size(Iapp));
for k = 1:size(time,2)
    time(:,k) = halfcycData(k).time(:);
end

% Compute average SOC and lithiation of porous electrode as a function of time.
[Q, theta0, theta100] = getCellParams(...
    simData.cellModel,'const.Q pos.theta0 pos.theta100','Output','list');
Z0 = soc0Pct/100; % initial SOC [fractional]
Zavg = zeros(size(time)); % average SOC vs time [fractional]
for k = 1:size(Zavg,2)
    Zavg(:,k) = Z0(k) - cumtrapz(time(:,k),Iapp(:,k))/Q/3600;
end
ThetaAvg = theta0 + Zavg*(theta100-theta0);

% Get Vcell, Iapp over common lithiation.
thetaAvg = linspace(min(ThetaAvg,[],'all'),max(ThetaAvg,[],'all'),1000).';
IappNorm = zeros(length(thetaAvg),size(Iapp,2));
VcellNorm = zeros(length(thetaAvg),size(Vcell,2));
for k = 1:size(Iapp,2)
    IappNorm(:,k) = interp1(ThetaAvg(:,k),Iapp(:,k),thetaAvg,'linear','extrap');
    VcellNorm(:,k) = interp1(ThetaAvg(:,k),Vcell(:,k),thetaAvg,'linear','extrap');
end

% Estimate OCP and resistance at each lithiation point using least-squares
% linear regression.
UocpEst = zeros(size(thetaAvg));
RcellEst = zeros(size(thetaAvg));
for k = 1:size(IappNorm,1)
    % Collect all current magnitudes applies and cell voltages measured 
    % at this lithiation point into column vectors.
    iapp = IappNorm(k,:).';
    vcell = VcellNorm(k,:).';
    % Compute least-squares solution to the linear system: 
    % vcell(theta) = Uocp(theta)- Rcell(theta)*iapp(theta)
    Ieqls0 = abs(iapp)<Q/200;  % logical indicies to zero current
    if sum(~Ieqls0)<2
        % Not enough nonzero current data-points to compute Uocp and Rcell
        % using LS; approximate Uocp only.
        UocpEst(k) = mean(vcell(Ieqls0));
        RcellEst(k) = NaN;
        continue;
    end
    H = [ones(size(iapp)) -iapp];  % measurement matrix
    x = H\vcell;
    UocpEst(k) = x(1);
    RcellEst(k) = x(2);
end

% Compute true OCP and Rcell for comparison.
LLPM = convertCellModel(simData.cellModel,'LLPM');
Rc = getCellParams(LLPM,'const.Rc');
ocpData = MSMR(LLPM.function.pos).ocp('theta',thetaAvg,'TdegC',simData.TdegC);
UocpTrue = ocpData.Uocp;
UocpRMSE = rms(UocpEst-UocpTrue);
resData = getPerturbationResistance(LLPM,thetaAvg);
RcellTrue = resData.Rtotal + Rc; % correct for tab resistance!
RcellIsNaN = isnan(RcellEst);
RcellRMSE = rms(RcellEst(~RcellIsNaN)-RcellTrue(~RcellIsNaN));

figure;
plot(thetaAvg,VcellNorm); hold on;
plot(thetaAvg,UocpTrue-IappNorm.*RcellTrue,':');
thesisFormat;

figure; 
plot(thetaAvg,UocpEst); hold on;
plot(thetaAvg,UocpTrue,':');
thesisFormat;

figure;
plot(thetaAvg,RcellEst); hold on;
plot(thetaAvg,RcellTrue,':');
thesisFormat;

return;

% -- Plotting -------------------------------------------------------------

% Plotting constants.
titlePrefix = sprintf( ...
    'Half-Cycle Discharge (%d\\rightarrow%d%% SOC)', ...
    simData.soc0Pct,simData.socfPct);
plotdir = fullfile('plots',simName);
if ~isfolder(plotdir)
    mkdir(plotdir);
end

% Plot Iapp(t).
labels1 = arrayfun( ...
    @(I)sprintf('%.1fC avg',I),simData.IavgC,'UniformOutput',false);
figure; colororder(spring(size(Iapp,2)));
plot(time/3600,Iapp/Q);
legend(labels1{:},'NumColumns',1,'Location','best');
title([titlePrefix ': i_{app}(t)']);
xlabel('Time, t [hr]');
ylabel('i_{app} [C-rate]');
thesisFormat([0.2 0.1 0.1 0.2]);

% Plot Vcell(t).
labels1 = arrayfun( ...
    @(I)sprintf('%.1fC avg FOM',I),simData.IavgC,'UniformOutput',false);
labels2 = arrayfun( ...
    @(err)sprintf('PM (%.3f%% RMSE)',err),rmsePct,'UniformOutput',false);
labels3 = arrayfun( ...
    @(I)sprintf('OCV'),simData.IavgC,'UniformOutput',false);
figure; colororder(spring(size(Vcell,2)));
plot(time/3600,Vcell); hold on;
plot(time/3600,VcellModel,'k--');
%plot(time/3600,Uocv,':');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','northeast');
%legend(labels1{:},labels2{:},labels3{:},'NumColumns',3,'Location','best');
title([titlePrefix ': v_{cell}(t)']);
xlabel('Time, t [hr]');
ylabel('v_{cell} [V]');
thesisFormat([0.2 0.1 0.1 0.2]);