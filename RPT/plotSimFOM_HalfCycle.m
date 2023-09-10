% plotSimFOM_HalfCycle.m
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
simName = 'cellSionGuess-P2DM-100pct-20pct';
load(fullfile('simdata','halfcyc',[simName '.mat']));

% Fetch time, iapp(t), vcell(t) vectors.
halfcycData = [simData.iappSeries];
Iapp = [halfcycData.iapp];
Vcell = [halfcycData.vcell];
time = zeros(size(Iapp));
for k = 1:size(time,2)
    time(:,k) = halfcycData(k).time(:);
end

% Compute average lithiation of porous electrode as a function
% of time (needed as an input to the perturbation model).
[Q, theta0, theta100] = getCellParams(...
    simData.cellModel,'const.Q pos.theta0 pos.theta100','Output','list');
Z0 = simData.soc0Pct/100; % initial SOC [fractional]
Zavg = zeros(size(time)); % average SOC vs time [fractional]
for k = 1:size(Zavg,2)
    Zavg(:,k) = Z0 - cumtrapz(time(:,k),Iapp(:,k))/Q/3600;
end
ThetaAvg = theta0 + Zavg*(theta100-theta0);

% Compute perturbation resistance and dynamic-equilibrium solution.
Rtotal = zeros(size(ThetaAvg));
Uocv = zeros(size(ThetaAvg));
clear parts;
for k = size(ThetaAvg,2):-1:1
    data = getPerturbationResistance( ...
        simData.cellModel,ThetaAvg(:,k),'TdegC',simData.TdegC);
    Rtotal(:,k) = data.Rtotal;
    parts(k) = data.parts;
    Uocv(:,k) = data.U;
end

% Compute Rctj over soc.
soc1 = linspace(5,100,1000);
thetaAvg1 = theta0 + (soc1/100)*(theta100-theta0);
data1 = getPerturbationResistance( ...
    simData.cellModel,thetaAvg1, ...
    'TdegC',simData.TdegC,'ComputeRctj',true);

% Compute vcell under perturbation model.
VcellModel = Uocv - Iapp.*Rtotal;

% Compute RMSE at each average C-rate.
rmsePct = sqrt(mean(((Vcell-VcellModel)./Vcell).^2,1))*100;

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
exportgraphics(gcf,fullfile(plotdir,'iapp.eps'));
exportgraphics(gcf,fullfile(plotdir,'iapp.png'));

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
exportgraphics(gcf,fullfile(plotdir,'vcell.eps'));
exportgraphics(gcf,fullfile(plotdir,'vcell.png'));

return;

% Plot Rtotal(t).
labels1 = arrayfun( ...
    @(I)sprintf('%.1fC avg',I),simData.IavgC,'UniformOutput',false);
figure; colororder(spring(size(Iapp,2)));
plot(time/3600,Rtotal); hold on;
legend(labels1{:},'NumColumns',1,'Location','best');
title([titlePrefix ': R_{total}(t)']);
xlabel('Time, t [hr]');
ylabel('R_{total} [\Omega]');
thesisFormat([0.2 0.1 0.1 0.2]);
exportgraphics(gcf,fullfile(plotdir,'Rtotal.eps'));
exportgraphics(gcf,fullfile(plotdir,'Rtotal.png'));

% Plot Rdiff(t).
labels1 = arrayfun( ...
    @(I)sprintf('%.1fC avg',I),simData.IavgC,'UniformOutput',false);
figure; colororder(spring(size(Iapp,2)));
plot(time/3600,[parts.Rdiff]);
legend(labels1{:},'NumColumns',1,'Location','best');
title([titlePrefix ': R_{diff}(t)']);
xlabel('Time, t [hr]');
ylabel('R_{diff} [\Omega]');
thesisFormat([0.2 0.1 0.1 0.2]);
exportgraphics(gcf,fullfile(plotdir,'Rdiff.eps'));
exportgraphics(gcf,fullfile(plotdir,'Rdiff.png'));

% Plot Rctp(t).
labels1 = arrayfun( ...
    @(I)sprintf('%.1fC avg',I),simData.IavgC,'UniformOutput',false);
figure; colororder(spring(size(Iapp,2)));
plot(time/3600,[parts.Rct_p]);
legend(labels1{:},'NumColumns',1,'Location','best');
title([titlePrefix ': R_{ct}^p(t)']);
xlabel('Time, t [hr]');
ylabel('R_{ct}^p [\Omega]');
thesisFormat([0.2 0.1 0.1 0.2]);
exportgraphics(gcf,fullfile(plotdir,'Rctp.eps'));
exportgraphics(gcf,fullfile(plotdir,'Rctp.png'));

% Plot Rctpj(t) vs soc.
J = simData.cellModel.metadata.section.pos.ocp.J;
labels1 = arrayfun(@(j)sprintf('R_{ct,%d}^p',j),1:J,'UniformOutput',false);
figure; colororder(spring(J));
semilogy(soc1,data1.parts.Rctj_p,':'); hold on;
semilogy(soc1,data1.parts.Rct_p,'r-');
yline(data.parts.Rct_n,'b-','LineWidth',3);
ylim([ ...
    min([data1.parts.Rct_p(:);data.parts.R0]) ...
    max([data1.parts.Rct_p(:);data.parts.R0]) ...
]);
set(gca,'xdir','reverse');
xlabel('Cell SOC [%]');
ylabel('Resistance Component [\Omega]');
title('Charge-Transfer Resistances vs. SOC');
legend(labels1{:},'R_{ct}^p','R_{ct}^n','NumColumns',3,'Location','north');
thesisFormat([0.2 0.1 0.1 0.2]);
exportgraphics(gcf,fullfile(plotdir,'Rctjp.eps'));
exportgraphics(gcf,fullfile(plotdir,'Rctjp.png'));