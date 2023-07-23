% fitPR.m
%
% Regress pulse models to resistance surfaces derived from laboratory
% measurments.
%
% -- Changelog --
% 08.17.2022 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath('..');
TB.addpaths;
addpath(TB.const.OCPROOT);

% Load resistance surfaces.
pr.load('SionFreshTMP-5tau-100000us','surfaces');
surface = surfaces.getTest(25); % select temperature in degC

% Load OCP model.
ocp.load('FinalFit-SionFresh_0C01','fit','varname','ocpmodels');
msmr = ocpmodels.getTest(40).model;  % 40degC is most accurate model

% Define optimization bounds.
lb.pos.k0 = 1e-6;       ub.pos.k0 = 1e3;
lb.pos.alpha = 0.4;     ub.pos.alpha = 0.6;
lb.pos.Rf = 1e-3;       ub.pos.Rf = 10;
lb.pos.Rdl = 1e-3;      ub.pos.Rdl = 10;
lb.pos.kappa = 1e-2;    ub.pos.kappa = 100;
lb.pos.sigma = 1;       ub.pos.sigma = 1e4;
lb.neg.k0 = 1e-3;       ub.neg.k0 = 10;
lb.neg.alpha = 0.4;     ub.neg.alpha = 0.6;
lb.neg.Rdl = 1e-3;      ub.neg.Rdl = 10;
% Cell tab resistance lumps Rf(n), 1/kappa(s), and 1/kappa(d).
lb.const.Rc = 0.01;     ub.const.Rc = 100;

% Perform the regression.
[result, modelspec, trajectory] = pr.fitR0( ...
    surface,msmr,lb,ub, ...
    'particleCount',500,'swarmIterations',500, ...
    'fminconIterations',10000, ...
    'trackTrajectory',true);
save( ...
    fullfile(TB.const.SDPROOT,'fit',sprintf('%s.mat',surface.name)), ...
    'result','modelspec','lb','ub','trajectory','msmr','surface');