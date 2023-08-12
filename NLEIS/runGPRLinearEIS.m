% runGPRLinearEIS.m
%
% Produce estimates of solid diffusivity and charge-transfer resistance
% using Gaussian Process Regression (GPR).
%
% 2023.07.29 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath('..');
TB.addpaths;
addpath(genpath('gpml'));
filename = 'EIS-16degC26degC-Ds=linear-k0=linear';
fitData = load(fullfile('labfitdata',[filename '.mat']));
plotdir = fullfile('plots','linEIS-GPR');
if ~isfolder(plotdir)
    mkdir(plotdir);
end

% Solid diffusivity -------------------------------------------------------
theta = fitData.values.pos.DsTheta;
Ds = fitData.values.pos.DsLinear;

% Define input and target data-points.
xtr = theta(:);      % training inputs
ytr = log10(Ds(:));  % training targets
xte = linspace(0,1,1000)';  % test inputs

% Define the prior as a GP with mean and covariance functions.
meanFn = {@meanConst};  % constant mean
covFn = {@covSEiso};    % squared exponential (SE) covariance
c = mean(ytr);  % mean value
sf = 1;         % signal std (std on f without noise)
ell = 0.05;     % length scale in input units

% Define the measurement liklihood.
likFn = {@likGauss};    % Gaussian liklihood distribution
sn = 0.2;       % measurement noise std

% Define initial hyperparameter vectors.
hyp0.mean = mean(ytr);
hyp0.cov = [log(ell); log(sf)];
hyp0.lik = log(sn);

% Minimize marginal liklihood over hyperparameters.
hyp = minimize(hyp0,@gp,-1000,@infGaussLik,meanFn,covFn,likFn,xtr,ytr);

% Compute prdictions.
[mu0,Sigma0] = gp(hyp0,@infGaussLik,meanFn,covFn,likFn,xtr,ytr,xte);
[mu,Sigma] = gp(hyp,@infGaussLik,meanFn,covFn,likFn,xtr,ytr,xte);

% Separate task: optimize extended-MSMR (e-MSMR) model.
params.Dsref = fastopt.param('logscale',true);
params.mD = fastopt.param;
spec = fastopt.modelspec(params);
lb.Dsref = 1e-9;     ub.Dsref = 1;
lb.mD = 1;           ub.mD = 5;
init.Dsref = 1e-5;
init.mD = 3;
lb = fastopt.pack(lb,spec);
ub = fastopt.pack(ub,spec);
init = fastopt.pack(init,spec);
electrode = MSMR(fitData.values.pos);
costFn = msmrDiffusivityRegressionCostFnFactory( ...
    theta,Ds,electrode,spec,fitData.arg.TrefdegC ...
);
vect = fmincon(costFn,init,[],[],[],[],lb,ub);
msmrDiffusionModel = fastopt.unpack(vect,spec);
postOpt_eMSMR_DsData = electrode.Ds(msmrDiffusionModel, ...
    'thetamin',0.01,'thetamax',0.99,'TdegC',fitData.arg.TrefdegC);

% Compute prediction of eMSMR model.
fitData = load('labfitdata\EIS-16degC26degC-Ds=msmr-k0=linear.mat');
eMSMR_DsData = electrode.Ds(fitData.values.pos, ...
    'thetamin',0.01,'thetamax',0.99,'TdegC',fitData.arg.TrefdegC);

% Compute prediction of 7-spline model.
fitData = load('labfitdata\EIS-16degC26degC-Ds=spline-k0=linear.mat');
splineDsData = electrode.Ds(fitData.values.pos, ...
    'thetamin',0.01,'thetamax',0.99,'TdegC',fitData.arg.TrefdegC);
thetaSpline = fitData.values.pos.DsTheta;
DsSpline = fitData.values.pos.DsSpline;

figure;
semilogy(xte,10.^mu0); hold on;
fill([xte;flipud(xte)],10.^[mu0+3*sqrt(Sigma0);flipud(mu0-3*sqrt(Sigma0))],...
       'k','EdgeColor','k','FaceAlpha',0.1,'EdgeAlpha',0.3);
semilogy(theta,Ds,'o');
set(gca,'xdir','reverse');
xlabel('$x$ in $\mathrm{Li}_x\mathrm{Ni}_y\mathrm{Mn}_z\mathrm{Co}_{1-y-z}\mathrm{O}_2$', ...
    'Interpreter','latex');
ylabel('Extrinsic Diffusivity, D_s [s^{-1}]');
title('Solid Diffusivity: Initial GPR Estimate');
legend('Estimate','3\sigma Bounds','From Linear EIS','Location','best');
thesisFormat;
print('-depsc',fullfile(plotdir,'Ds-initial'));
print('-dpng',fullfile(plotdir,'Ds-initial'));

figure;
semilogy(xte,10.^mu); hold on;
fill([xte;flipud(xte)],10.^[mu+3*sqrt(Sigma);flipud(mu-3*sqrt(Sigma))],...
       'k','EdgeColor','k','FaceAlpha',0.1,'EdgeAlpha',0.3);
semilogy(theta,Ds,'o');
set(gca,'xdir','reverse');
xlabel('$x$ in $\mathrm{Li}_x\mathrm{Ni}_y\mathrm{Mn}_z\mathrm{Co}_{1-y-z}\mathrm{O}_2$', ...
    'Interpreter','latex');
ylabel('Extrinsic Diffusivity, D_s [s^{-1}]');
title('Solid Diffusivity: Optimized GPR Estimate');
legend('Estimate','3\sigma Bounds','From Linear EIS','Location','best');
thesisFormat;
print('-depsc',fullfile(plotdir,'Ds-optimized'));
print('-dpng',fullfile(plotdir,'Ds-optimized'));

figure;
semilogy(postOpt_eMSMR_DsData.theta,postOpt_eMSMR_DsData.Ds); hold on;
semilogy(theta,Ds,'o');
set(gca,'xdir','reverse');
xlabel('$x$ in $\mathrm{Li}_x\mathrm{Ni}_y\mathrm{Mn}_z\mathrm{Co}_{1-y-z}\mathrm{O}_2$', ...
    'Interpreter','latex');
ylabel('Extrinsic Diffusivity, D_s [s^{-1}]');
title('Solid Diffusivity: Optimized eMSMR Estimate');
legend('Estimate','From Linear EIS','Location','best');
thesisFormat;
print('-depsc',fullfile(plotdir,'Ds-eMSMR-optimized'));
print('-dpng',fullfile(plotdir,'Ds-eMSMR-optimized'));

figure;
semilogy(xte,10.^mu); hold on;
semilogy(eMSMR_DsData.theta,postOpt_eMSMR_DsData.Ds);
semilogy(splineDsData.theta,splineDsData.Ds);
set(gca,'xdir','reverse');
xlabel('$x$ in $\mathrm{Li}_x\mathrm{Ni}_y\mathrm{Mn}_z\mathrm{Co}_{1-y-z}\mathrm{O}_2$', ...
    'Interpreter','latex');
ylabel('Extrinsic Diffusivity, D_s [s^{-1}]');
title('Solid Diffusivity: Comparing Estimates');
legend('GPR','eMSMR','7-spline','Location','best');
thesisFormat;
print('-depsc',fullfile(plotdir,'Ds-compare'));
print('-dpng',fullfile(plotdir,'Ds-compare'));

% Exchange-current --------------------------------------------------------
theta = fitData.values.pos.k0Theta;
k0 = fitData.values.pos.k0Linear;

% Define input and target data-points.
xtr = theta(:);            % training inputs
ytr = log10(k0(:));        % training targets
xte = linspace(0,1,1000)'; % test inputs

% Define the prior as a GP with mean and covariance functions.
meanFn = {@meanConst};  % constant mean
covFn = {@covSEiso};    % squared exponential (SE) covariance
c = mean(ytr);  % mean value
sf = 0.5;         % signal std (std on f without noise)
ell = 0.02;     % length scale in input units

% Define the measurement liklihood.
likFn = {@likGauss};   % Gaussian liklihood distribution
sn = 0.1;              % measurement noise std

% Define initial hyperparameter vectors.
hyp0.mean = mean(ytr);
hyp0.cov = [log(ell); log(sf)];
hyp0.lik = log(sn);

% Minimize marginal liklihood over hyperparameters.
hyp = minimize(hyp0,@gp,-1000,@infGaussLik,meanFn,covFn,likFn,xtr,ytr);

% Compute prdictions.
[mu0,Sigma0] = gp(hyp0,@infGaussLik,meanFn,covFn,likFn,xtr,ytr,xte);
[mu,Sigma] = gp(hyp,@infGaussLik,meanFn,covFn,likFn,xtr,ytr,xte);

figure;
semilogy(xte,10.^mu0); hold on;
fill([xte;flipud(xte)],10.^[mu0+3*sqrt(Sigma0);flipud(mu0-3*sqrt(Sigma0))],...
       'k','EdgeColor','k','FaceAlpha',0.1,'EdgeAlpha',0.3);
semilogy(theta,k0,'o');
set(gca,'xdir','reverse');
xlabel('$x$ in $\mathrm{Li}_x\mathrm{Ni}_y\mathrm{Mn}_z\mathrm{Co}_{1-y-z}\mathrm{O}_2$', ...
    'Interpreter','latex');
ylabel('Exchange Current, i_{0}^{p} [A]');
title('i_{o}^{p}: Initial GPR Estimate');
legend('Estimate','3\sigma Bounds','From Linear EIS','Location','best');
thesisFormat;
print('-depsc',fullfile(plotdir,'i0p-initial'));
print('-dpng',fullfile(plotdir,'i0p-initial'));

figure;
semilogy(xte,10.^mu); hold on;
fill([xte;flipud(xte)],10.^[mu+3*sqrt(Sigma);flipud(mu-3*sqrt(Sigma))],...
       'k','EdgeColor','k','FaceAlpha',0.1,'EdgeAlpha',0.3);
semilogy(theta,k0,'o');
set(gca,'xdir','reverse');
xlabel('$x$ in $\mathrm{Li}_x\mathrm{Ni}_y\mathrm{Mn}_z\mathrm{Co}_{1-y-z}\mathrm{O}_2$', ...
    'Interpreter','latex');
ylabel('Exchange Current , i_{0}^{p} [A]');
title('i_{0}^{p}: Optimized GPR Estimate');
legend('Estimate','3\sigma Bounds','From Linear EIS','Location','best');
thesisFormat;
print('-depsc',fullfile(plotdir,'i0p-optimized'));
print('-dpng',fullfile(plotdir,'i0p-optimized'));

function fn = msmrDiffusivityRegressionCostFnFactory( ...
    theta,DsLab,electrode,modelspec,TdegC ...
)
    ocpData = electrode.ocp('theta',theta,'TdegC',TdegC);
    fn = @cost;

    function J = cost(vect)
        model = fastopt.unpack(vect,modelspec);
        dsData = electrode.DsCachedOCP(model,ocpData);
        J = sum((log(DsLab(:))-log(dsData.Ds(:))).^2);
    end % cost()
end % costFnFactory()