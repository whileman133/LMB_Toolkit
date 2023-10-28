% runGPRLinearEIS.m
%
% Produce estimates of solid diffusivity and charge-transfer resistance
% using Gaussian Process Regression (GPR).
%
% 2023.08.27 | Employ a number of intermediate solns in GPR | Wes H
% 2023.07.29 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath('..');
TB.addpaths;
addpath(genpath('gpml'));
filename = '202309_EIS-16degC26degC-Ds=linear-k0=linear';
fitData = load(fullfile('labfitdata',[filename '.mat']));
plotdir = fullfile('plots','linEIS-GPR');
if ~isfolder(plotdir)
    mkdir(plotdir);
end

% Constants.
NumIntermediateSolns = 20;

% Solid diffusivity -------------------------------------------------------
% Fetch observations of Ds over lithiation theta. Include the "best"
% solution and also similar (but not exactly equal) solutions. We ensure
% the solutions differ by computing the normalized distance between each
% similar solution and the "best" solution.
pos = [fitData.topSoln.Ds(1:NumIntermediateSolns).pos];
theta = [pos.DsTheta];
Ds = [pos.DsLinear];
theta = theta(:);
Ds = Ds(:);

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
hyp0.mean = c;
hyp0.cov = [log(ell); log(sf)];
hyp0.lik = log(sn);

% Minimize marginal liklihood over hyperparameters.
hyp = minimize(hyp0,@gp,-1000,@infGaussLik,meanFn,covFn,likFn,xtr,ytr);

% Compute prdictions.
[mu0,Sigma0] = gp(hyp0,@infGaussLik,meanFn,covFn,likFn,xtr,ytr,xte);
[mu,Sigma] = gp(hyp,@infGaussLik,meanFn,covFn,likFn,xtr,ytr,xte);
est.log10Ds = mu;
est.sigma = Sigma;
est.theta = xte;
save(fullfile('gpr',[filename '_DsEstimate.mat']),'est');

% Separate task: see if we can fit MSMR model to this.
params.Dsref = fastopt.param('logscale',true);
params.zeta = fastopt.param();
spec = fastopt.modelspec(params);
lb.Dsref = 1e-9;  ub.Dsref = 100;
lb.zeta  = 0;   ub.zeta  = 1;
init.Dsref = 1e-5;
init.zeta = 0.5;
lb = fastopt.pack(lb,spec);
ub = fastopt.pack(ub,spec);
init = fastopt.pack(init,spec);
electrode = MSMR(fitData.values.pos);
costFn = msmrDiffusivityRegressionCostFnFactory( ...
    fitData.values.pos.DsTheta,fitData.values.pos.DsLinear, ...
    electrode,spec,fitData.arg.TrefdegC ...
);
vect = fmincon(costFn,init,[],[],[],[],lb,ub);
msmrDiffusionModel = fastopt.unpack(vect,spec);
zeta = msmrDiffusionModel.zeta;
Dsref1 = msmrDiffusionModel;
postOpt_MSMR_DsData = electrode.Ds(msmrDiffusionModel, ...
     'thetamin',0.01,'thetamax',0.99,'TdegC',fitData.arg.TrefdegC);

% Separate task: see if we can fit (extended) MSMR model to this.
clear params spec lb ub init electrode;
params.Dsref = fastopt.param('logscale',true);
params.theta0 = fastopt.param('fix',fitData.values.pos.theta0);
params.theta100 = fastopt.param('fix',fitData.values.pos.theta100);
params.mD = fastopt.param();
spec = fastopt.modelspec(params);
lb.Dsref = 1e-9;  ub.Dsref = 100;
lb.theta0 = 0.9;  ub.theta0 = 1;
lb.theta100 = 0;  ub.theta100 = 0.4;
lb.mD = 1;        ub.mD = 10;
init.Dsref = 1e-5;
init.theta0 = fitData.values.pos.theta0;
init.theta100 = fitData.values.pos.theta100;
init.mD = 2.5;
lb = fastopt.pack(lb,spec);
ub = fastopt.pack(ub,spec);
init = fastopt.pack(init,spec);
electrode = MSMR(fitData.values.pos);
costFn = msmrDiffusivityRegressionCostFnFactoryVariableTheta( ...
    fitData.values.pos.DsTheta,fitData.values.pos.DsLinear, ...
    electrode,spec,fitData.arg.TrefdegC ...
);
vect = fmincon(costFn,init,[],[],[],[],lb,ub);
msmrDiffusionModel = fastopt.unpack(vect,spec);
postOpt_eMSMR_DsData = electrode.Ds(msmrDiffusionModel, ...
     'thetamin',0.01,'thetamax',0.99,'TdegC',fitData.arg.TrefdegC);
soc = (fitData.values.pos.DsTheta-fitData.values.pos.theta0)/(fitData.values.pos.theta100-fitData.values.pos.theta0);
theta_eMSMR = msmrDiffusionModel.theta0 + soc.*(msmrDiffusionModel.theta100-msmrDiffusionModel.theta0);

% Separate task: optimize extended-MSMR (e-MSMR) model.
% params.Dsref = fastopt.param('logscale',true);
% params.mD = fastopt.param;
% spec = fastopt.modelspec(params);
% lb.Dsref = 1e-9;     ub.Dsref = 1;
% lb.mD = 1;           ub.mD = 5;
% init.Dsref = 1e-5;
% init.mD = 3;
% lb = fastopt.pack(lb,spec);
% ub = fastopt.pack(ub,spec);
% init = fastopt.pack(init,spec);
% electrode = MSMR(fitData.values.pos);
% costFn = msmrDiffusivityRegressionCostFnFactory( ...
%     theta,Ds,electrode,spec,fitData.arg.TrefdegC ...
% );
% vect = fmincon(costFn,init,[],[],[],[],lb,ub);
% msmrDiffusionModel = fastopt.unpack(vect,spec);
% postOpt_eMSMR_DsData = electrode.Ds(msmrDiffusionModel, ...
%     'thetamin',0.01,'thetamax',0.99,'TdegC',fitData.arg.TrefdegC);

% Compute prediction of eMSMR model.
% fitData = load('labfitdata\EIS-16degC26degC-Ds=msmr-k0=linear.mat');
% eMSMR_DsData = electrode.Ds(fitData.values.pos, ...
%     'thetamin',0.01,'thetamax',0.99,'TdegC',fitData.arg.TrefdegC);

% Compute prediction of 7-spline model.
% fitData = load('labfitdata\EIS-16degC26degC-Ds=spline-k0=linear.mat');
% splineDsData = electrode.Ds(fitData.values.pos, ...
%     'thetamin',0.01,'thetamax',0.99,'TdegC',fitData.arg.TrefdegC);
% thetaSpline = fitData.values.pos.DsTheta;
% DsSpline = fitData.values.pos.DsSpline;

% Plot initial GPR estimate.
% figure;
% semilogy(xte,10.^mu0); hold on;
% fill([xte;flipud(xte)],10.^[mu0+3*sqrt(Sigma0);flipud(mu0-3*sqrt(Sigma0))],...
%        'k','EdgeColor','k','FaceAlpha',0.1,'EdgeAlpha',0.3);
% semilogy(theta,Ds,'o');
% set(gca,'xdir','reverse');
% xlabel('$x$ in $\mathrm{Li}_x\mathrm{Ni}_y\mathrm{Mn}_z\mathrm{Co}_{1-y-z}\mathrm{O}_2$', ...
%     'Interpreter','latex');
% ylabel('Extrinsic Diffusivity, D_s [s^{-1}]');
% title('Solid Diffusivity: Initial GPR Estimate');
% legend('Estimate','3\sigma Bounds','From Linear EIS','Location','best');
% thesisFormat;
% print('-depsc',fullfile(plotdir,'Ds-initial'));
% print('-dpng',fullfile(plotdir,'Ds-initial'));

% Plot optimized GPR estimate.
figure;
semilogy(xte,10.^mu,'k:'); hold on;
fill([xte;flipud(xte)],10.^[mu+3*sqrt(Sigma);flipud(mu-3*sqrt(Sigma))],...
       'k','EdgeColor','k','FaceAlpha',0.1,'EdgeAlpha',0.3);
semilogy(theta,Ds,'ro');
semilogy(xte,10.^mu,'k:');
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
semilogy(xte,10.^mu,'k:'); hold on;
fill([xte;flipud(xte)],10.^[mu+3*sqrt(Sigma);flipud(mu-3*sqrt(Sigma))],...
       'k','EdgeColor','k','FaceAlpha',0.1,'EdgeAlpha',0.3);
semilogy(postOpt_MSMR_DsData.theta,postOpt_MSMR_DsData.Ds);
set(gca,'xdir','reverse');
xlabel('$x$ in $\mathrm{Li}_x\mathrm{Ni}_y\mathrm{Mn}_z\mathrm{Co}_{1-y-z}\mathrm{O}_2$', ...
    'Interpreter','latex');
ylabel('Extrinsic Diffusivity, D_s [s^{-1}]');
title('Solid Diffusivity: Optimized MSMR Estimate');
legend('GPR Estimate','3\sigma Bounds','MSMR Estimate','Location','northwest');
thesisFormat;
print('-depsc',fullfile(plotdir,'Ds-MSMR-optimized'));
print('-dpng',fullfile(plotdir,'Ds-MSMR-optimized'));

figure;
semilogy(xte,10.^mu,'k:'); hold on;
fill([xte;flipud(xte)],10.^[mu+3*sqrt(Sigma);flipud(mu-3*sqrt(Sigma))],...
       'k','EdgeColor','k','FaceAlpha',0.1,'EdgeAlpha',0.3);
semilogy(postOpt_eMSMR_DsData.theta,postOpt_eMSMR_DsData.Ds);
set(gca,'xdir','reverse');
xlabel('$x$ in $\mathrm{Li}_x\mathrm{Ni}_y\mathrm{Mn}_z\mathrm{Co}_{1-y-z}\mathrm{O}_2$', ...
    'Interpreter','latex');
ylabel('Extrinsic Diffusivity, D_s [s^{-1}]');
title('Solid Diffusivity: Optimized e-MSMR Estimate');
legend('GPR Estimate','3\sigma Bounds','MSMR Estimate','Location','northwest');
thesisFormat;
print('-depsc',fullfile(plotdir,'Ds-eMSMR-optimized'));
print('-dpng',fullfile(plotdir,'Ds-eMSMR-optimized'));

return;

% figure;
% semilogy(postOpt_eMSMR_DsData.theta,postOpt_eMSMR_DsData.Ds); hold on;
% semilogy(theta,Ds,'o');
% set(gca,'xdir','reverse');
% xlabel('$x$ in $\mathrm{Li}_x\mathrm{Ni}_y\mathrm{Mn}_z\mathrm{Co}_{1-y-z}\mathrm{O}_2$', ...
%     'Interpreter','latex');
% ylabel('Extrinsic Diffusivity, D_s [s^{-1}]');
% title('Solid Diffusivity: Optimized eMSMR Estimate');
% legend('Estimate','From Linear EIS','Location','best');
% thesisFormat;
% print('-depsc',fullfile(plotdir,'Ds-eMSMR-optimized'));
% print('-dpng',fullfile(plotdir,'Ds-eMSMR-optimized'));

% figure;
% semilogy(xte,10.^mu); hold on;
% semilogy(eMSMR_DsData.theta,postOpt_eMSMR_DsData.Ds);
% semilogy(splineDsData.theta,splineDsData.Ds);
% set(gca,'xdir','reverse');
% xlabel('$x$ in $\mathrm{Li}_x\mathrm{Ni}_y\mathrm{Mn}_z\mathrm{Co}_{1-y-z}\mathrm{O}_2$', ...
%     'Interpreter','latex');
% ylabel('Extrinsic Diffusivity, D_s [s^{-1}]');
% title('Solid Diffusivity: Comparing Estimates');
% legend('GPR','eMSMR','7-spline','Location','best');
% thesisFormat;
% print('-depsc',fullfile(plotdir,'Ds-compare'));
% print('-dpng',fullfile(plotdir,'Ds-compare'));

% Exchange-current --------------------------------------------------------
pos = [fitData.topSoln.k0(1:NumIntermediateSolns).pos];
theta = [pos.k0Theta];
k0 = [pos.k0Linear];
theta = theta(:);
k0 = k0(:);
%theta = fitData.values.pos.k0Theta;
%k0 = fitData.values.pos.k0Linear;

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

function fn = msmrDiffusivityRegressionCostFnFactoryVariableTheta( ...
    theta,DsLab,electrode,modelspec,TdegC ...
)
    soc = (theta-electrode.zmax)/(electrode.zmin-electrode.zmax);
    fn = @cost;

    function J = cost(vect)
        model = fastopt.unpack(vect,modelspec);
        t = model.theta0 + soc*(model.theta100-model.theta0);
        model.lm = (model.theta100-model.theta0)/(electrode.zmin-electrode.zmax);
        ocpData = electrode.ocp('theta',t,'TdegC',TdegC);
        dsData = electrode.DsCachedOCP(model,ocpData);
        J = sum((log(DsLab(:))-log(dsData.Ds(:))).^2);
    end % cost()
end % costFnFactory()