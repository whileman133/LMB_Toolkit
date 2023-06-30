function regressionData = fitLinearEIS(spectra,ocptest,ocpmodel)
%FITLINEAREIS

% Constants
tplus0 = 0.4;    % Li+ transference number, guess for estimating psi
R = TB.const.R;  % molar gas constant [J/mol/K]
F = TB.const.F;  % Faraday's constant [C/mol]

% Build model -------------------------------------------------------------
% Known parameters.
params.const.Q = fastopt.param('fix',ocptest.QAh);
params.pos.theta0 = fastopt.param('fix',ocpmodel.zmax);
params.pos.theta100 = fastopt.param('fix',ocpmodel.zmin);
params.pos.X = fastopt.param('fix',ocpmodel.X);
params.pos.U0 = fastopt.param('fix',ocpmodel.U0);
params.pos.omega = fastopt.param('fix',ocpmodel.omega);

% Fixed parameters.
% Rdl(p|n) are low for the Sion cells, so fix Rdl(p|n)=0.
params.pos.Rdl = fastopt.param('fix',0);
params.neg.Rdl = fastopt.param('fix',0);
params.neg.Rf = fastopt.param('fix',0);  % Lump into tab resistance.
% Setting wDL=0 results in ideal CPEs for the double-layer capacitors;
% okay since we will not be converting to the time domain inside of the
% optimization routine.
params.pos.wDL = fastopt.param('fix',0);
params.neg.wDL = fastopt.param('fix',0);
% Symmetry factors for the positive electrode are not very identifiable,
% we'll asume alpha=0.5 for all galleries.
params.pos.alpha = fastopt.param('fix',0.5*ones(ocpmodel.J,1));
% Symmetry factor for the negative electrode is not identifiable due to
% linear nature of the small-signal impedance model.
params.neg.alpha = fastopt.param('fix',0.5);

% Electrolyte parameters.
% psi,W,tauW are not separately identifiable, so we fix psi=R/(F*(1-t+0)), 
% the value predicted by the Einstien relationship.
params.const.psi = fastopt.param('fix',R/F/(1-tplus0));
params.const.W = fastopt.param;
params.pos.tauW = fastopt.param;
params.pos.kappa = fastopt.param;
params.eff.tauW = fastopt.param;
params.eff.kappa = fastopt.param;

% Porous electrode parameters.
params.pos.k0 = fastopt.param('fix',zeros(ocpmodel.J,1),'fixmask',k0mask);
params.pos.Dsref = fastopt.param('logscale',true);
params.pos.nF = fastopt.param;
params.pos.sigma = fastopt.param;
params.pos.Cdl = fastopt.param;
params.pos.nDL = fastopt.param;
params.pos.Rf = fastopt.param;

% Lithium-metal electrode parameters.
params.neg.k0 = fastopt.param;
params.neg.Cdl = fastopt.param;
params.neg.nDL = fastopt.param;

% Cell package parameters.
params.pkg.R0 = fastopt.param;  % Tab resistance.
params.pkg.L0 = fastopt.param;  % Package/cable inductace.

modelspec = fastopt.modelspec(params);


end