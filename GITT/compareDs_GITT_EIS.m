% compareDs_GITT_EIS.m

clear; close all; clc;
addpath('..');
TB.addpaths;

eisFile = '202309_EIS-16degC26degC-Ds=linear-k0=linear_DsEstimate.mat';
gittFiles = {
    'GITT-1.mat'
    'GITT-4.mat'
    'GITT-5.mat'
};

% Load data files.
eisData = load(fullfile(TB.const.NLEISROOT,'gpr',eisFile));
clear gittData;
for i = length(gittFiles):-1:1
    gittData(i) = load(fullfile('processed',gittFiles{i}));
end % for

% Scale relative Ds estimates from GITT to match EIS as best as possible 
% using least squares regression.
DsEIS = 10.^eisData.est.log10Ds(:);
ksquared = zeros(size(gittData));
for i = 1:length(gittData)
    est = gittData(i).est;
    DsGITT = 10.^est.log10Ds(:);
    C = [log10(DsGITT) ones(length(DsGITT),1)];
    d = log10(DsEIS);
    Aeq = [1 0];
    beq = 1;
    coeff = lsqlin(C,d,[],[],Aeq,beq);
    ksquared(i) = 10.^coeff(2);
end % for

% Compare estimates.
figure;
line = semilogy(eisData.est.theta,10.^eisData.est.log10Ds,'DisplayName','Linear EIS'); hold on;
colorEIS = line.Color;
colorsGITT = zeros(length(gittData),3);
for i = 1:length(gittData)
    name = split(gittFiles{i},'.');
    name = name{1};
    ksq = ksquared(i);
    est = gittData(i).est;
    line = semilogy(est.theta,(10.^est.log10Ds)*ksq,'DisplayName',sprintf('%s (k^2=%.2e m/s)',name,ksq));
    colorsGITT(i,:) = line.Color;
end
fill([eisData.est.theta;flipud(eisData.est.theta)], ...
    10.^[eisData.est.log10Ds+3*sqrt(eisData.est.sigma);flipud(eisData.est.log10Ds-3*sqrt(eisData.est.sigma))],...
   colorEIS,'EdgeColor',colorEIS,'FaceAlpha',0.1,'EdgeAlpha',0.3,'DisplayName','3\sigma Bounds');
for i = 1:length(gittData)
    name = split(gittFiles{i},'.');
    name = name{1};
    ksq = ksquared(i);
    est = gittData(i).est;
    color = colorsGITT(i,:);
    fill([est.theta;flipud(est.theta)], ...
        10.^[est.log10Ds+3*sqrt(est.Sigma);flipud(est.log10Ds-3*sqrt(est.Sigma))]*ksq,...
       color,'EdgeColor',color,'FaceAlpha',0.1,'EdgeAlpha',0.3,'DisplayName','3\sigma Bounds');
end
set(gca,'xdir','reverse');
xlabel('$x$ in $\mathrm{Li}_x\mathrm{Ni}_y\mathrm{Mn}_z\mathrm{Co}_{1-y-z}\mathrm{O}_2$', ...
    'Interpreter','latex');
ylabel('Diffusion Coefficient, $D_s$ [$s^{-1}$]','Interpreter','latex');
legend('Location','best','NumColumns',2);
title('Comparing Diffusion Coeff. from EIS and GITT');
thesisFormat;
print(fullfile('plots','compare','all'),'-depsc');
print(fullfile('plots','compare','all'),'-dpng');