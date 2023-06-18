% makeModelLMO.m
%
% Make a lumped-parameter model describing a full LMB cell.
% LMO porous electrode.
%
% -- Changelog --
% 01.22.2023 | Created | Wesley Hileman <whileman@uccs.edu>
%
% -- References --
% We use parameter values from this literature:
% [1] Shanshan Xu et al 2019 J. Electrochem. Soc. 166 A3456
% [2] Daniel R. Baker and Mark W. Verbrugge 2021 J. Electrochem. Soc. 168 050526
% [3] Dongliang Lu et al 2022 J. Electrochem. Soc. 169 080504

clear; close all; clc;

% Define parameters of P2D model.
p.const.A = 64e-4;      % cell cross-sectional area [m^2]
p.const.F = 96485.3321; % Faraday constant [C/mol]
p.const.R = 8.31446262; % molar gas constant [J/mol·K]
p.const.dlnfdlnc = 3;   % from Lu et al. [3]
p.const.t0plus = 0.4;   % from Baker and Verbrugge [2] [unitless]
p.const.kappa = 0.8;    % @(T,ce0) from Baker and Verbrugge [2] [1/Ω·m]
p.const.De = 2.58e-9;   % @(T,ce0) from Baker and Verbrugge [2] [m^2/s]
p.const.T = 298;        % from Baker and Verbrugge [2] [K]
p.const.ce0 = 1000;     % from Baker and Verbrugge [2] [mol/m^3]
p.const.brug = 1.5;     % from Lu et al. [3] [unitless]
p.pos.L = 200e-6;       % from Baker and Verbrugge [2] [m]
p.pos.epsilons = 0.44;  % from Baker and Verbrugge [2] [m]
p.pos.epsilone = 0.37;  % from Baker and Verbrugge [2] [m]
p.pos.Rs = 4e-6;        % from Baker and Verbrugge [2] [m]
p.pos.as = 3*p.pos.epsilons/p.pos.Rs;
p.pos.sigma = 25.9;     % from Baker and Verbrugge [2] [1/Ω·m]
p.pos.Dsref = 1.33e-14; % from Baker and Verbrugge [2] [m^2/s]
p.pos.csmax = 20000;    % from Baker and Verbrugge [2] [mol/m^3]
p.pos.nF = 1;           % [unitless]
p.pos.nDL = 1;          % [unitless]
p.pos.Rf = 0.02;        % [Ω·m^2]
p.pos.Rdl = 0;          % [Ω·m^2]
p.pos.Cdl = 0.22;       % [F/m^2]
p.pos.U0    = [4.16756 4.02477];    % from Baker and Verbrugge [2] [V]
p.pos.X     = [1-0.60331 0.60331];  % from Baker and Verbrugge [2] [unitless]
p.pos.omega = [1.12446 1.71031];    % from Baker and Verbrugge [2] [unitless]
p.pos.alpha = [0.2 0.71307];        % from Baker and Verbrugge [2] [unitless]
p.pos.k0    = [145.6 132.1]/10;     % from Baker and Verbrugge [2] [A/m^2]
p.pos.theta0 = 0.99;
p.pos.theta100 = 0.01;
p.sep.L = 25e-6;        % from Baker and Verbrugge [2] [m]
p.sep.epsilon = 0.41;   % from Baker and Verbrugge [2] [unitless]
p.sep.nE = 1;           % [unitless]
p.DL.L = 5e-6;          % [m]
p.DL.epsilon = 0.21;    % [unitless]
p.DL.nE = 1;            % [unitless]
p.neg.k0 = 3.5e-8;      % from Xu et al [1] [m/s]
p.neg.cs = 76900;       % from Xu et al [1] [mol/m^3]
p.neg.alpha = 0.7;      % from Xu et al [1] [unitless]
p.neg.gamma = 1;        % [unitless]
p.neg.nDL = 1;          % [unitless]
p.neg.Rf = 0.02;        % [Ω·m^2]
p.neg.Rdl = 0;          % [Ω·m^2]
p.neg.Cdl = 0.198;      % [F/m^2]

% Convert P2D model to lumped form.
l.const.F = p.const.F;
l.const.R = p.const.R;
l.const.Q = p.pos.epsilons*p.const.A*p.pos.L*p.const.F*p.pos.csmax*abs(p.pos.theta100-p.pos.theta0)/3600; % [Ah]
l.const.psi = p.const.F*p.const.De*p.const.ce0/p.const.kappa/(1-p.const.t0plus)/p.const.T; % [V/K]
l.const.kD = 2*p.const.R*(p.const.t0plus-1)*(1+p.const.dlnfdlnc)/p.const.F; % [V/K]
l.pos.sigma = p.const.A*p.pos.sigma*(p.pos.epsilons)^p.const.brug/p.pos.L; % [1/Ω]
l.pos.kappa = p.const.A*p.const.kappa*(p.pos.epsilone)^p.const.brug/p.pos.L; % [1/Ω]
l.pos.qe = p.const.F*p.pos.epsilone*p.const.ce0*p.const.A*p.pos.L/(1-p.const.t0plus)/3600; % [Ah]
l.pos.Dsref = p.pos.Dsref/p.pos.Rs^2; % [1/s]
l.pos.nF = p.pos.nF; % [unitless]
l.pos.nDL = p.pos.nDL; % [unitless]
l.pos.Rf = p.pos.Rf/p.pos.as/p.const.A/p.pos.L; % [Ω]
l.pos.Rdl = p.pos.Rdl/p.pos.as/p.const.A/p.pos.L; % [Ω]
l.pos.Cdl = p.pos.Cdl*p.pos.as*p.const.A*p.pos.L; % [F]
l.pos.U0 = p.pos.U0; % [V]
l.pos.X = p.pos.X; % [unitless]
l.pos.omega = p.pos.omega; % [unitless]
l.pos.alpha = p.pos.alpha; % [unitless]
l.pos.theta0 = p.pos.theta0;
l.pos.theta100 = p.pos.theta100;
l.pos.k0 = p.const.A*p.pos.L*p.pos.as*p.pos.k0; % [A]
l.sep.kappa = p.const.A*p.const.kappa*(p.sep.epsilon)^p.const.brug/p.sep.L; % [1/Ω]
l.sep.qe = p.const.F*p.sep.epsilon*p.const.ce0*p.const.A*p.sep.L/(1-p.const.t0plus)/3600; % [Ah]
l.sep.nE = p.sep.nE;
l.DL.kappa = p.const.A*p.const.kappa*(p.DL.epsilon)^p.const.brug/p.DL.L; % [1/Ω]
l.DL.qe = p.const.F*p.DL.epsilon*p.const.ce0*p.const.A*p.DL.L/(1-p.const.t0plus)/3600; % [Ah]
l.DL.nE = p.DL.nE;
l.neg.k0 = p.neg.gamma*p.const.A*p.const.F*p.neg.k0*(p.neg.cs)^p.neg.alpha*p.const.ce0^(1-p.neg.alpha); % [A]
l.neg.alpha = p.neg.alpha; % [unitless]
l.neg.nDL = p.neg.nDL; % [unitless]
l.neg.Rf = p.neg.Rf/p.neg.gamma/p.const.A; % [Ω]
l.neg.Rdl = p.neg.Rdl/p.neg.gamma/p.const.A; % [Ω]
l.neg.Cdl = p.neg.Cdl*p.neg.gamma*p.const.A; % [F]

% Save model.
lumped = l; p2d = p;
save('cellLMO.mat','lumped','p2d');