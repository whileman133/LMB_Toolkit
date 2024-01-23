% demoScaleOCP.m
%
% Show how thetamin and thetamax scale the OCP curve of a porous electrode.
%
%
% -- Changelog --
% 2024.01.21 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath('..');
TB.addpaths;

% Fetch OCP over absolute composition for NMC622 electrode.
electrode = MSMR.NMC622();
ocp = electrode.ocp('npoints',100);
% Compute relative composition of electrode.
thetaTilde = (ocp.theta-electrode.zmin)/(electrode.zmax-electrode.zmin);

figure;
plot(thetaTilde,ocp.Uocp,':'); hold on;
plot(ocp.theta,ocp.Uocp,'-');
xline(electrode.zmin,'k--');
xline(electrode.zmax,'k--');
text(electrode.zmin+0.005,min(ocp.Uocp)-0.02,'zmin', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment','left', ...
    'FontSize',16);
text(electrode.zmax-0.0051,min(ocp.Uocp)-0.02,'zmax', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment','right', ...
    'FontSize',16);
xlabel('theta');
ylabel('Uocp');
title('ScaleOCP');
legend('Uocp1','Uocp2','FontSize',24);
thesisFormat;
print(fullfile('plots','scaleOCP'),'-depsc');
print(fullfile('plots','scaleOCP'),'-dpng');