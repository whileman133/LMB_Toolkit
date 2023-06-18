% makeAgedModelLMO.m
%
% Make lumped-parameter models describing a full LMB cell over the course
% of cycle aging. LMO porous electrode.
%
% -- Changelog --
% 2023.06.06 | Created | Wesley Hileman <whileman@uccs.edu>
%
% -- References --
% We use parameter values from this literature:
% [1] Shanshan Xu et al 2019 J. Electrochem. Soc. 166 A3456
% [2] Daniel R. Baker and Mark W. Verbrugge 2021 J. Electrochem. Soc. 168 050526
% [3] Dongliang Lu et al 2022 J. Electrochem. Soc. 169 080504

clear; close all; clc;

% Normalized variable that describes state-of-age of the cell.
ageVect = linspace(0,1,10);

% Baseline model.
p.const.Rc = 0;         % tab/contact resistance [Ω]
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

clear std lumped;
for k = length(ageVect):-1:1
    age = ageVect(k);  % age variable
    aged = p;          % aged model

    % Electrolyte depletion.
    aged.const.ce0 = aged.const.ce0*(1-age*0.5);      % ce0  ..  0.5ce0
    aged.const.kappa = aged.const.kappa*(1-age*0.7);  % k    ..  0.3k
    aged.const.De = aged.const.De*(1-age*0.7);        % De   ..  0.3De

    % Dendrites and dead-lithium accumulation at negative electrode.
    aged.neg.gamma = aged.neg.gamma*(1+age*2);    % gam  ..  3gam
    aged.neg.Rf = aged.neg.Rf*(1+age*4);          % Rf   ..  5Rf
    aged.DL.L = aged.DL.L*(1+age*9);              % L    ..  10L
    aged.DL.epsilon = aged.DL.epsilon*(1-age/2);  % eps  ..  0.5eps

    % Particle cracking at positive electrode.
    aged.pos.Rs = aged.pos.Rs*(1-age*0.9);        % Rs   ..  0.1Rs   
    aged.pos.as = 3*aged.pos.epsilons/p.pos.Rs;
    aged.pos.Rf = aged.pos.Rf*(1+age*4);          % Rf   ..  5Rf

    % Phase change in porous-electrode lattice.
    aged.pos.Dsref = aged.pos.Dsref*(1-age*0.8);  % Ds   ..  0.2Ds
    aged.pos.U0 = [
        aged.pos.U0(1)*(1+age*0.10)   % 0%  ..  +10% drift
        aged.pos.U0(2)*(1+age*0.07)   % 0%  ..  +7% drift
    ];
    aged.pos.omega = aged.pos.omega*(1+age*0.1);  % 0%  ..  +10% drift

    % Change in reaction symmetry.
    aged.neg.alpha = aged.neg.alpha*(1-age*0.3);  % 0%  ..  -30% drift
    aged.pos.alpha = [
        aged.pos.alpha(1)*(1+age*0.3)   % 0%  ..  +30% drift
        aged.pos.alpha(2)*(1-age*0.3)   % 0%  ..  -30% drift
    ];

    lump = lumpCellModel(aged);
    lump.function = functionify(lump);
    lump.function.pos.soc = eval( ...
        sprintf('@(x)(%g+x*%g)', ...
            lump.pos.theta0,lump.pos.theta100-lump.pos.theta0) ...
    );

    std(k) = aged;
    lumped(k) = lump;
end

% Save model.
save('cellLMO_AgeSeries.mat','lumped','std','ageVect');


function fcn = functionify(model)
    paramNames = fieldnames(model);
    for k = 1:length(paramNames)
        paramName = paramNames{k};
        value = model.(paramName);
        if isstruct(value)
            fcn.(paramName) = functionify(value);
        elseif isscalar(value)
            % Scalar value.
            fcnString = sprintf('@(x,T)(%g)',value);
            fcn.(paramName) = eval(fcnString);
        else
            % Vector value.
            fcnString = ['@(x,T)([' sprintf('%g;',value(:)) '])'];
            fcn.(paramName) = eval(fcnString);
        end
    end
end