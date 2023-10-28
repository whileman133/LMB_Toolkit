function data = gprBasic(xtr,ytr,xte,hyp0)
%GPR Perform Gaussian process regression (GPR) with constant mean and 
%  squared exponential covariance.
%
% data = gprBasic(xtr, ytr, xte, sf, ell, sn)
%  xtr  .. training inputs
%  ytr  .. training targets
%  xte  .. test inputs
%  hyp0.sf  .. initial signal std (std on f without noise) in output units
%  hyp0.ell .. initial squared-exponential length scale in input units
%  hyp0.sn  .. initial measurement noise std in output units

% Force column vectors.
xtr = xtr(:);
ytr = ytr(:);
xte = xte(:);

% Define the prior as a GP with mean and covariance functions.
meanFn = {@meanConst};  % constant mean
covFn = {@covSEiso};    % squared exponential (SE) covariance
likFn = {@likGauss};    % Gaussian liklihood distribution
c = mean(ytr);  % mean value

% Define initial hyperparameter vectors.
hy0.mean = c;
hy0.cov = [log(hyp0.ell); log(hyp0.sf)];
hy0.lik = log(hyp0.sn);

% Minimize marginal liklihood over hyperparameters.
hy = minimize(hy0,@gp,-1000,@infGaussLik,meanFn,covFn,likFn,xtr,ytr);

% Compute prdictions.
[mu0,Sigma0] = gp(hy0,@infGaussLik,meanFn,covFn,likFn,xtr,ytr,xte);
[mu,Sigma] = gp(hy,@infGaussLik,meanFn,covFn,likFn,xtr,ytr,xte);

data.initial.mu = mu0;
data.initial.Sigma = Sigma0;
data.optimized.mu = mu;
data.optimized.Sigma = Sigma;

end