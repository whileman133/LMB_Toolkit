% plotSimFOM_MultipleHalfCycle.m
%
% Plot simulated half-cycle discharge, compare to Baker and Verbrugge
% perturbation approximation.
%
% -- Changelog --
% 2023.12.20 | Use processSimHalfCycle.m | Wesley Hileman
% 2023.06.01 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath(fullfile('..','..'));
TB.addpaths();

% Load simulation data.
simName = 'cellLMO-P2DM-99pct-1pct';
load(fullfile('simdata',[simName '.mat']));
plotdir = fullfile('plots',simName);
cellName = simData.cellModel.metadata.cell.name;

[theta0, theta100] = getCellParams(...
    simData.cellModel,'pos.theta0 pos.theta100','Output','list');
hcycData = processSimHalfCycle(simData.iappSeries);
UocpEst = hcycData.ocp.Uocp;
UocpSocAvgPct = hcycData.ocp.socAvgPct;
UocpThetaAvg = theta0 + (UocpSocAvgPct/100)*(theta100-theta0);
RdcEst = hcycData.res.Rdc;
RdcSocAvgPct = hcycData.res.socAvgPct;
RdcThetaAvg = theta0 + (RdcSocAvgPct/100)*(theta100-theta0);
VcellNorm = hcycData.VcellNorm;
IappNorm = hcycData.IappNorm;

% Compute true OCP and Rcell for comparison.
LLPM = convertCellModel(simData.cellModel,'LLPM');
ocpData = MSMR(LLPM.function.pos).Ds( ...
    LLPM.function.pos,'theta',UocpThetaAvg,'TdegC',simData.TdegC);
UocpTrue = ocpData.Uocp;
UocpRMSE = rms(UocpEst-UocpTrue);
resData = getDcResistance(LLPM,RdcThetaAvg);
RcellTrue = resData.Rdc;
RcellRMSE = rms(RdcEst(~isnan(RdcEst))-RcellTrue(~isnan(RdcEst)));

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
plot(UocpThetaAvg,VcellNorm); hold on;
plot(UocpThetaAvg,VcellHat,':');
set(gca,'xdir','reverse');
xlabel('Fractional Lithium Composition, $\theta_\mathrm{s}$', ...
    'Interpreter','latex');
ylabel('Cell voltage, v_{cell} [V]');
title(sprintf('Half-Cycle Dis/charge (%s)',cellName));
l = legend(lab1{:},lab2{:},'Location','','NumColumns',2);
thesisFormat;
l.FontSize = 10;
l.Position(1) = 0.2;
l.Position(2) = 0.32;
print(fullfile(plotdir,'hc-dis'),'-depsc');
print(fullfile(plotdir,'hc-dis'),'-dpng');

figure;
% subplot(211);
lab = arrayfun(@(I)sprintf('Iavg=%+.2fC',I),simData.IavgC,'UniformOutput',false);
plot(UocpThetaAvg,(VcellHat-VcellTrue)*1000); hold on;
set(gca,'xdir','reverse');
xlabel('Fractional Lithium Composition, $\theta_\mathrm{s}$', ...
    'Interpreter','latex');
ylabel('Prediction error, $\hat{v}_{cell}-v_{cell}$ [mV]','Interpreter','latex');
title(sprintf('Prediction Error of dc Model'));
legend(lab,'Location','best','NumColumns',1);
% subplot(212); 
% yyaxis left;
% plot(UocpThetaAvg,ocpData.d2Uocp);
% set(gca,'xdir','reverse');
% xlabel('Fractional Lithium Composition, $\theta_\mathrm{s}$', ...
%     'Interpreter','latex');
% ylabel('OCP curvature, U_{ocp}'''' [V]');
% title('OCP Curvature and Solid Diffusivity');
% yyaxis right;
% plot(UocpThetaAvg,ocpData.Ds,':');
% ylabel('Solid-diffusion coefficient D_s [s^{-1}]');
thesisFormat('PlotBoxPaddingInches',[0 0 0 0]);
print(fullfile(plotdir,'vcell-error'),'-depsc');
print(fullfile(plotdir,'vcell-error'),'-dpng');

figure; 
plot(UocpThetaAvg,UocpTrue); hold on;
plot(UocpThetaAvg,UocpEst,':');
set(gca,'xdir','reverse');
xlabel('Fractional Lithium Composition, $\theta_\mathrm{s}$', ...
    'Interpreter','latex');
ylabel('U_{ocp} [V vs. Li/Li+]');
title(sprintf('OCP Estimate (%.3fmV RMSE)',UocpRMSE*1000));
legend('True OCP','Estimate','Location','best');
thesisFormat;
print(fullfile(plotdir,'OCPEst'),'-depsc');
print(fullfile(plotdir,'OCPEst'),'-dpng');

figure;
plot(RdcThetaAvg,RcellTrue); hold on;
plot(RdcThetaAvg,RdcEst,':');
set(gca,'xdir','reverse');
xlabel('Fractional Lithium Composition, $\theta_\mathrm{s}$', ...
    'Interpreter','latex');
ylabel('R_{dc} [\Omega]');
title(sprintf('DC Resistance Estimate (%.3fm\\Omega RMSE)',RcellRMSE*1000));
legend('True DC Resistance','Estimate','Location','best');
thesisFormat;
print(fullfile(plotdir,'RdcEst'),'-depsc');
print(fullfile(plotdir,'RdcEst'),'-dpng');
