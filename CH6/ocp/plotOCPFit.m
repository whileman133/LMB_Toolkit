% plotOCPFit.m

clear; close all; clc;
addpath(fullfile('..','..'));
TB.addpaths;

% Constants ---------------------------------------------------------------

% Name of the study to load.
STUDYNAME = 'Fit-SionFresh_0C01';  % C/100 quasi-equilibrium cycling

% Temperature select.
TDEGC = 40;


% Load data ---------------------------------------------------------------

ocp.load(STUDYNAME,'fit','studydir',fullfile('.','labdata','fit'));
fit = fitstudy.getTest(TDEGC);
p.U0 = fit.model.Uj0;
p.X = fit.model.Xj;
p.omega = fit.model.Wj;
p.thetamin = fit.model.zmin;
p.thetamax = fit.model.zmax;
model = MSMR(p);

% Move laboratory measurements over absolute composition.
lab = fit.ocpest;
lababs.Z = p.thetamin + lab.Z*(p.thetamax-p.thetamin);
lababs.V = lab.V;
lababs.refV = lab.refV;
lababs.refZ = p.thetamin + lab.refZ*(p.thetamax-p.thetamin);
lababs.dzdvRefZ = lab.dzdvRefZ*(p.thetamax-p.thetamin);
lababs.dzdvRefV = lab.dzdvRefV*(p.thetamax-p.thetamin);

% Plotting ----------------------------------------------------------------

ocpm = model.ocp('theta',lababs.Z);
ocpRMSE = rms(ocpm.Uocp(:)-lababs.V(:));

figure;
plot(ocpm.theta,ocpm.Uocp); hold on;
plot(lababs.Z,lababs.V,':','Color',[0 0.8 0]);
ylim([min(ocpm.Uocp) max(ocpm.Uocp)]);
legend('Model------------','Lab','FontSize',20,'Location','Southwest');
title('MSMRreg');
xlabel('theta');
ylabel('Uocp');
thesisFormat;
print(fullfile('plots','OCP-fit'),'-depsc');
print(fullfile('plots','OCP-fit'),'-dpng');

ocpm = model.ocp('theta',lababs.refZ);
diffcapRMSE = rms(-1./ocpm.dUocp(:)-lababs.dzdvRefZ(:));

figure;
plot(ocpm.theta,abs(1./ocpm.dUocp)); hold on;
plot(lababs.refZ,lababs.dzdvRefZ,':','Color',[0 0.8 0]);
xlim([0 1]);
legend('Model------------','Lab','FontSize',20,'Location','Northwest');
title('MSMRreg');
xlabel('theta');
ylabel('diffcap');
thesisFormat;
print(fullfile('plots','diffcap-Z-fit'),'-depsc');
print(fullfile('plots','diffcap-Z-fit'),'-dpng');

ocpm = model.ocp('voltage',lababs.refV);

figure;
plot(abs(1./ocpm.dUocp),ocpm.Uocp); hold on;
plot(lababs.dzdvRefV,lababs.refV,':','Color',[0 0.8 0]);
ylim([min(ocpm.Uocp) max(ocpm.Uocp)]);
legend('Model------------','Lab','FontSize',20,'Location','Southeast');
title('MSMRreg');
xlabel('diffcap');
ylabel('Uocp');
thesisFormat;
print(fullfile('plots','diffcap-V-fit'),'-depsc');
print(fullfile('plots','diffcap-V-fit'),'-dpng');