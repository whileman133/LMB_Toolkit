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
simName = 'cellSionGuess-P2DM-99pct-1pct';
load(fullfile('simdata','multhalfcyc',[simName '.mat']));
plotdir = fullfile('plots','multhalfcycle',simName);

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
ocpData = MSMR(LLPM.function.pos).Ds( ...
    LLPM.function.pos,'theta',thetaAvg,'TdegC',simData.TdegC);
UocpTrue = ocpData.Uocp;
UocpRMSE = rms(UocpEst-UocpTrue);
resData = getDcResistance(LLPM,thetaAvg);
RcellTrue = resData.Rdc;
RcellIsNaN = isnan(RcellEst);
RcellRMSE = rms(RcellEst(~RcellIsNaN)-RcellTrue(~RcellIsNaN));

% Compute voltage-prediction RMSE.
VcellTrue = VcellNorm;
VcellHat = UocpTrue-IappNorm.*RcellTrue;
VcellRMSE = rms(VcellTrue-VcellHat);
disp(simData.IavgC);
disp(VcellRMSE*1000);

% -- Plotting -------------------------------------------------------------
if ~isfolder(plotdir)
    mkdir(plotdir);
end

figure;
lab1 = arrayfun(@(I)sprintf('COMSOL'),simData.IavgC,'UniformOutput',false);
lab2 = arrayfun(@(I)sprintf('DC Model @ %+.2fC',I),simData.IavgC,'UniformOutput',false);
plot(thetaAvg,VcellNorm); hold on;
plot(thetaAvg,UocpTrue-IappNorm.*RcellTrue,':');
set(gca,'xdir','reverse');
xlabel('$x$ in $\mathrm{Li}_x\mathrm{Ni}_y\mathrm{Mn}_z\mathrm{Co}_{1-y-z}\mathrm{O}_2$', ...
    'Interpreter','latex');
ylabel('Cell voltage, v_{cell} [V]');
title('Half-Cycle Dis/charge');
l = legend(lab1{:},lab2{:},'Location','best','NumColumns',2);
thesisFormat;
l.FontSize = 10;
print(fullfile(plotdir,'hc-dis'),'-depsc');
print(fullfile(plotdir,'hc-dis'),'-dpng');

figure;
subplot(211);
lab = arrayfun(@(I)sprintf('Iavg=%+.2fC',I),simData.IavgC,'UniformOutput',false);
plot(thetaAvg,(VcellHat-VcellTrue)*1000); hold on;
set(gca,'xdir','reverse');
xlabel('$x$ in $\mathrm{Li}_x\mathrm{Ni}_y\mathrm{Mn}_z\mathrm{Co}_{1-y-z}\mathrm{O}_2$', ...
    'Interpreter','latex');
ylabel('Prediction error, $\hat{v}_{cell}-v_{cell}$ [mV]','Interpreter','latex');
title('Prediction Error of dc Model');
legend(lab,'Location','best','NumColumns',2);
subplot(212); 
yyaxis left;
plot(thetaAvg,ocpData.d2Uocp);
set(gca,'xdir','reverse');
xlabel('$x$ in $\mathrm{Li}_x\mathrm{Ni}_y\mathrm{Mn}_z\mathrm{Co}_{1-y-z}\mathrm{O}_2$', ...
    'Interpreter','latex');
ylabel('OCP curvature, U_{ocp}'''' [V]');
title('OCP Curvature and Solid Diffusivity');
yyaxis right;
plot(thetaAvg,ocpData.Ds,':');
ylabel('Solid-diffusion coefficient D_s [s^{-1}]');
thesisFormat('PlotBoxPaddingInches',[0 0 0.5 0]);
print(fullfile(plotdir,'d2Uocp-Ds'),'-depsc');
print(fullfile(plotdir,'d2Uocp-Ds'),'-dpng');
return;

figure; 
plot(thetaAvg,UocpTrue); hold on;
plot(thetaAvg,UocpEst,':');
set(gca,'xdir','reverse');
xlabel('$x$ in $\mathrm{Li}_x\mathrm{Ni}_y\mathrm{Mn}_z\mathrm{Co}_{1-y-z}\mathrm{O}_2$', ...
    'Interpreter','latex');
ylabel('U_{ocp} [V vs. Li/Li+]');
title(sprintf('OCP Estimate (%.3fmV RMSE)',UocpRMSE*1000));
legend('True OCP','Estimate','Location','best');
thesisFormat;
print(fullfile(plotdir,'OCPEst'),'-depsc');
print(fullfile(plotdir,'OCPEst'),'-dpng');

figure;
plot(thetaAvg,RcellTrue); hold on;
plot(thetaAvg,RcellEst,':');
set(gca,'xdir','reverse');
xlabel('$x$ in $\mathrm{Li}_x\mathrm{Ni}_y\mathrm{Mn}_z\mathrm{Co}_{1-y-z}\mathrm{O}_2$', ...
    'Interpreter','latex');
ylabel('R_{dc} [\Omega]');
title(sprintf('DC Resistance Estimate (%.3fm\\Omega RMSE)',RcellRMSE*1000));
legend('True DC Resistance','Estimate','Location','best');
thesisFormat;
print(fullfile(plotdir,'RdcEst'),'-depsc');
print(fullfile(plotdir,'RdcEst'),'-dpng');
