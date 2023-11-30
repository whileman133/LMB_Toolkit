% runFitRPTSim.m
% 
% Regress cell model to RPT data collected in simulation to verify the RPT
% protocol.
%
% -- Changelog --
% 2023.11.14 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath('..');
TB.addpaths;

simData = load(fullfile('simdata','fullrpt','cellLMO_AgeArray.mat'));
meas = simData.ageSeries(end,end);
hcyc = meas.halfCycle;
eis = meas.eis;

ageData = load(fullfile('agemodel','cellLMO_AgeArray.mat'));
modelTrue = convertCellModel(ageData.ageArray(end,end).model,'WRM');
pTrue = getCellParams(modelTrue);

% Constants.
tplus0 = 0.4;    % Li+ transference number, guess for estimating psi
R = TB.const.R;  % molar gas constant [J/mol/K]
F = TB.const.F;  % Faraday's constant [C/mol]
% Theta vector for linear LUT for pos.k0
thetaLin = linspace(pTrue.pos.theta100,pTrue.pos.theta0,1000);

% Build model -------------------------------------------------------------
% Known parameters.
params.const.Q = fastopt.param('fix',pTrue.const.Q);
params.pos.theta0 = fastopt.param('fix',pTrue.pos.theta0);
params.pos.theta100 = fastopt.param('fix',pTrue.pos.theta100);
params.pos.X = fastopt.param('fix',pTrue.pos.X);
params.pos.U0 = fastopt.param('fix',pTrue.pos.U0);
params.pos.omega = fastopt.param('fix',pTrue.pos.omega);

% Fixed parameters.
% Rdl(p|n) are low for the Sion cells, so fix Rdl(p|n)=0.
params.pos.Rdl = fastopt.param('fix',0);
params.neg.Rdl = fastopt.param('fix',0);
params.neg.Rf = fastopt.param('fix',0);  % Lump into tab resistance.
% Symmetry factors for the positive electrode are not very identifiable,
% so we won't try to estimate them.
params.pos.alpha = fastopt.param('fix',pTrue.pos.alpha);
% Symmetry factor for the negative electrode is not identifiable due to
% linear nature of the small-signal impedance model.
params.neg.alpha = fastopt.param('fix',pTrue.neg.alpha);
% Solid conductance does not influence impedance significatly.
params.pos.sigma = fastopt.param('fix',pTrue.pos.sigma);
% psi,W,tauW are not separately identifiable, so we fix psi=R/(F*(1-t+0)), 
% the value predicted by the Einstien relationship.
params.const.psi = fastopt.param('fix',R/F/(1-tplus0));
% Fix lithiation of k0 interpolation points.
params.pos.k0Theta = fastopt.param('fix',thetaLin);

% Electrolyte parameters.
params.const.W = fastopt.param('logscale',true);
params.pos.tauW = fastopt.param('logscale',true);
params.pos.kappa = fastopt.param('logscale',true);
params.eff.tauW = fastopt.param('logscale',true);
params.eff.kappa = fastopt.param('logscale',true);

% Porous electrode parameters.
params.pos.Dsref = fastopt.param('logscale',true);
params.pos.mD = fastopt.param;
if strcmpi(arg.KineticsModel,'linear')
    params.pos.k0Linear = fastopt.param('len',length(thetaLin),'logscale',true,'tempfcn','Eact');
else
    params.pos.k0Spline = fastopt.param('len',ocpmodel.J,'logscale',true,'tempfcn','Eact');
end
params.pos.nF = fastopt.param;
params.pos.tauF = fastopt.param('logscale',true);
params.pos.Cdl = fastopt.param;
params.pos.nDL = fastopt.param;
params.pos.tauDL = fastopt.param('logscale',true);
params.pos.Rf = fastopt.param;

% Lithium-metal electrode parameters.
params.neg.k0 = fastopt.param('logscale',true,'tempfcn','Eact');
params.neg.Cdl = fastopt.param;
params.neg.nDL = fastopt.param;
params.neg.tauDL = fastopt.param('logscale',true);

% Cell package parameters.
params.pkg.R0 = fastopt.param('tempfcn','lut');  % Tab resistance.
params.pkg.L0 = fastopt.param('tempfcn','lut');  % Package/cable inductace.

modelspec = fastopt.modelspec(params, ...
    'tempsdegC',TdegC,'TrefdegC',arg.TrefdegC);