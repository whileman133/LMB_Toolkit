% plotUocpXPD.m
%
% Compares the NMC-811 OCP curve from Friedrich et al. to our laboratory
% measurements.
%
% 05.11.2022
% Wesley Hileman <whileman@uccs.edu>
% University of Colorado Colorado Springs

clear; close all; clc;
addpath(fullfile('..','..'));
TB.addpaths;

ocp.load('SionFresh_0C01','built');
load(fullfile('.','labdata','UocpXPD.mat'));
zXPD = x; UocpXPD = Uocp;

% Size of the voltage bins used for computing differential capacity.
VBINSIZE = 5e-3;

ocptest = builtstudy.getTest(40);
ocpest = ocp.DiagonalEstimate(ocptest,VBINSIZE,'maxZ',0.7,'minV',3.7).useAvg;
ocp = ocpest.asOCP();

% Absolute lithiation range corresponding to Uocp=[4.3,3.2]V range.
zmin = fzero(@(z)(interp1(x,UocpXPD,z,'linear','extrap')-4.3),[0 1]); % @ 4.3V
zmax = 1; % @ 3.2V (no data, so guess)

% When U=4.3V, theta=zmin. When U=3.2V, theta=zmax.
zocp = zmin + (zmax - zmin)*(ocp.ZU - ocp.zmin)/(ocp.zmax - ocp.zmin);
Uocp = ocp.U;

figure;
plot(zocp,Uocp); hold on;
plot(zXPD(1:120:end),UocpXPD(1:120:end),'.','Color',[0 0.8 0]);
xline(zmin,'k'); yline(4.3,'k'); xline(0.3651,'k--');
xlim([0 1]); ylim([3.2 4.6]);
annotation('textarrow',[0.21 0.21],[0.67 0.73],'String','VMAX','FontSize',20);
annotation('textarrow',[0.33 0.293],[0.3 0.3],'String','TMIN','FontSize',20);
annotation('textarrow',[0.49 0.45],[0.65 0.65],'String','BND','FontSize',20);
xlabel('theta');
ylabel('Potential');
legend('DiagonalAveraging-------------','Literature','FontSize',20);
title('XPD');
thesisFormat;
print(fullfile('plots','UocpXPD'),'-depsc');
print(fullfile('plots','UocpXPD'),'-dpng');

% Plot Baker-Verbrugge Diffusivity
% refZocp = zmin + (zmax - zmin)*(ocpest.refZ - ocp.zmin)/(ocp.zmax - ocp.zmin);
% dUocp = 1./ocpest.dzdvRefZ;
% Dnorm = MSMR.F*(refZocp).*(1-refZocp).*abs(dUocp)/MSMR.R/(ocpest.ocptest.temp+273.15);
% 
% figure;
% plot(refZocp,Dnorm); hold on;
% xlim([0 1]);
% xlabel('Absolute Lithiation, \theta_s');
% ylabel('Normalized Diffusivity, D_s/D_{s,ref}');
% title('Lithiation-Dependent Solid Diffusivity');
% thesisFormat([0.1 0.05 0.1 0.1]);
% todisk(gcf,fullfile('..','plots','BakerVerbruggeDsXPD'));
% 
% % Plot OCP relationship (scaling).
% figure;
% plot(ocp.ZU,ocp.U); hold on;
% plot(zocp,Uocp);
% xlim([0 1]); ylim([3.2 4.3]);
% xline(1-0.8221,'k--'); xline(1,'k--');
% annotation('textarrow',[0.35 0.275],[0.3 0.3],'String',' \theta_{min}','FontSize',20);
% annotation('textarrow',[0.89 0.955],[0.6 0.6],'String',' \theta_{max}','FontSize',20);
% xlabel('Absolute or Relative Lithiation');
% ylabel('Potential versus Li/Li^+ [V]');
% title('Scaling the OCP Curve');
% legend('Relative $\tilde{U}_{\mathrm{ocp}}(\tilde{\theta}_s)$ (Diagonal Avg.)', 'Scaled Absolute $U_{\mathrm{ocp}}(\theta_s)$', 'Interpreter','latex');
% thesisFormat([0.1 0.05 0.1 0.1]);
% todisk(gcf,fullfile('..','plots','UocpSPMScaling'));