% plotOCP.m
%
% Plot regressed MSMR OCP curve.
%
% -- Changelog --
% 07.23.2023 | Update for gen2 toolkit | Wesley Hileman
% 10.12.2022 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath('..');
TB.addpaths();
filename = 'FinalFit-SionFresh_0C01.mat';
TdegC = 40;
% Size of the voltage bins used for computing differential capacity.
% (for more information, see: help smoothdiff.m)
VBINSIZE = 10e-3;

load(fullfile(TB.const.OCPROOT,'labdata','fitstruct',filename));
load(fullfile(TB.const.OCPROOT,'labdata','fit',filename));

% Compute dU/d(theta) and d2U/d(theta)2 from regressed MSMR model.
test = study.tests(study.testTemperatures==TdegC);
electrode = MSMR(test.MSMR,'name',test.name);
ocpData = electrode.ocp('thetamin',0.001,'thetamax',0.999,'TdegC',TdegC);

% Compute dU/d(theta) and d2U/d(theta)2 from lab data.
zmax = electrode.zmax;
zmin = 0.18;  % due to augmented XPD data.
test = fitstudy.tests(fitstudy.testTemperatures==TdegC).ocpest.ocptest;
diagEst = ocp.DiagonalEstimate(test,VBINSIZE,'maxZ',0.7,'minV',3.7);
est = diagEst.useChg;
theta = zmin + est.Z*(zmax - zmin);
Uocp = est.V;
theta1 = zmin + est.refZ*(zmax - zmin);
dUocp = -1./est.dzdvRefZ./(zmax - zmin);
theta2 = theta1(1:end-1) + diff(theta1)/2;
d2Uocp = diff(dUocp)./diff(theta1);

% Plotting ----------------------------------------------------------------
ocp.MSMR.compare( ...
    3,5,2,test.temp, ...
    sprintf('Lab (%s %.0f\\circC)',test.name,test.temp),est, ...
    'Model Fit',ocp.MSMR(electrode.Xj,electrode.Uj0,electrode.Wj,zmin,zmax) ...
);
set(gcf,'Color',[239, 229, 195]/255);
set(gcf, 'InvertHardcopy', 'off');
thesisFormat;
print(fullfile('plots','Uocp-dUocp'),'-dpng');
print(fullfile('plots','Uocp-dUocp'),'-depsc');
return;

figure;
plot(ocpData.theta,ocpData.Uocp); hold on;
plot(theta,Uocp,':');
xlim([0.1 0.997]);
%set(gca,'xdir','reverse');
xlabel('$x$ in $\mathrm{Li}_x\mathrm{Ni}_y\mathrm{Mn}_z\mathrm{Co}_{1-y-z}\mathrm{O}_2$', ...
    'Interpreter','latex');
ylabel('Potential vs. $\mathrm{Li/Li}^+$, $U_\mathrm{ocp}$ [V]','Interpreter','latex');
legend('Lab','MSMR Model','Location','northeast');
title('OCP');
thesisFormat;
print(fullfile('plots','Uocp'),'-dpng');

figure;
plot(theta1,dUocp); hold on;
plot(ocpData.theta,ocpData.dUocp,':');
xlim([0.18 0.95]);
%set(gca,'xdir','reverse');
xlabel('$x$ in $\mathrm{Li}_x\mathrm{Ni}_y\mathrm{Mn}_z\mathrm{Co}_{1-y-z}\mathrm{O}_2$', ...
    'Interpreter','latex');
ylabel('OCP Slope, $\frac{\mathrm{d}U_\mathrm{ocp}}{\mathrm{d}x}$ [V]','Interpreter','latex');
legend('Lab','Model','Location','southeast');
title('OCP Slope');
thesisFormat;
print(fullfile('plots','dUocp'),'-dpng');

figure;
plot(theta2,d2Uocp); hold on;
plot(ocpData.theta,ocpData.d2Uocp,':');
xlim([0.18 0.95]);
%set(gca,'xdir','reverse');
xlabel('$x$ in $\mathrm{Li}_x\mathrm{Ni}_y\mathrm{Mn}_z\mathrm{Co}_{1-y-z}\mathrm{O}_2$', ...
    'Interpreter','latex');
ylabel('OCP Curvature, $\frac{\mathrm{d}^2U_\mathrm{ocp}}{\mathrm{d}x^2}$ [V]','Interpreter','latex');
legend('Lab','Model','Location','southwest');
title('OCP Curvature');
thesisFormat;
print(fullfile('plots','d2Uocp'),'-dpng');