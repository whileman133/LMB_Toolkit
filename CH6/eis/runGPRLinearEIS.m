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
filename = 'EIS-16degC26degC-Ds=linear-k0=linear';
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
init.mD = 1.5;
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
mD = msmrDiffusionModel.mD;
postOpt_eMSMR_DsData = electrode.Ds(msmrDiffusionModel, ...
     'thetamin',0.01,'thetamax',0.99,'TdegC',fitData.arg.TrefdegC);
soc = (fitData.values.pos.DsTheta-fitData.values.pos.theta0)/(fitData.values.pos.theta100-fitData.values.pos.theta0);
theta_eMSMR = msmrDiffusionModel.theta0 + soc.*(msmrDiffusionModel.theta100-msmrDiffusionModel.theta0);

% Plot optimized GPR estimate.
figure;
semilogy(xte,10.^mu,'k:'); hold on;
fill([xte;flipud(xte)],10.^[mu+3*sqrt(Sigma);flipud(mu-3*sqrt(Sigma))],...
       [0.9 0.9 0.9],'EdgeColor',[0.8 0.8 0.8],'FaceAlpha',1.0,'EdgeAlpha',1.0);
semilogy(theta,Ds,'ro');
semilogy(xte,10.^mu,'k:');
set(gca,'xdir','reverse');
xlabel('theta');
ylabel('Ds');
title('SolidDiff');
legend('GPR Mean','GPR 3\sigma Bounds','EIS','Location','best');
thesisFormat;
print('-depsc',fullfile(plotdir,'Ds-optimized'));
print('-dpng',fullfile(plotdir,'Ds-optimized'));

figure;
semilogy(xte,10.^mu,'k:'); hold on;
fill([xte;flipud(xte)],10.^[mu+3*sqrt(Sigma);flipud(mu-3*sqrt(Sigma))],...
       [0.9 0.9 0.9],'EdgeColor',[0.8 0.8 0.8],'FaceAlpha',1.0,'EdgeAlpha',1.0);
semilogy(postOpt_MSMR_DsData.theta,postOpt_MSMR_DsData.Ds);
semilogy(postOpt_eMSMR_DsData.theta,postOpt_eMSMR_DsData.Ds);
semilogy(xte,10.^mu,'k:'); hold on;
set(gca,'xdir','reverse');
xlabel('theta');
ylabel('Ds');
title('SolidDiff');
legend( ...
    'EIS (mean)','EIS (3\sigma Bounds)', ...
    'Baker-Verbrugge', ...
    sprintf('Stretched Baker-Verbrugge (m=%.1f)',mD), ...
    'Location','best');
thesisFormat;
print('-depsc',fullfile(plotdir,'Ds-MSMR-optimized'));
print('-dpng',fullfile(plotdir,'Ds-MSMR-optimized'));


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
semilogy(xte,10.^mu,'k:'); hold on;
fill([xte;flipud(xte)],10.^[mu+3*sqrt(Sigma);flipud(mu-3*sqrt(Sigma))],...
       [0.9 0.9 0.9],'EdgeColor',[0.8 0.8 0.8],'FaceAlpha',1.0,'EdgeAlpha',1.0);
semilogy(theta,k0,'ro');
semilogy(xte,10.^mu,'k:');
set(gca,'xdir','reverse');
xlabel('theta');
ylabel('i0');
title('ExchgCurrent');
legend('GPR Mean','GPR 3\sigma Bounds','EIS','Location','best');
thesisFormat;
print('-depsc',fullfile(plotdir,'i0p-optimized'));
print('-dpng',fullfile(plotdir,'i0p-optimized'));


% Utility Functions -------------------------------------------------------

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