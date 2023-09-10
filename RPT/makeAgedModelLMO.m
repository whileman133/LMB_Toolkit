% makeAgedModelLMO.m
%
% Make lumped-parameter models describing a full LMB cell over the course
% of cycle aging. LMO porous electrode.
%
% -- Changelog --
% 2023.09.07 | Update for gen2 toolkit | Wes H
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
model0 = loadCellModel('cellLMO-P2DM');

clear models;
for k = length(ageVect):-1:1
    age = ageVect(k);  % age variable
    aged = getCellParams(model0);  % initial structure of aged parameters

    % Electrolyte depletion.
    aged.const.ce0 = aged.const.ce0*(1-age*0.5);      % ce0  ..  0.5ce0
    aged.const.kappa = aged.const.kappa*(1-age*0.7);  % k    ..  0.3k
    aged.const.De = aged.const.De*(1-age*0.7);        % De   ..  0.3De

    % Dendrites and dead-lithium accumulation at negative electrode.
    aged.neg.gamma = aged.neg.gamma*(1+age*2);    % gam  ..  3gam
    aged.neg.Rf = aged.neg.Rf*(1+age*4);          % Rf   ..  5Rf
    aged.dll.L = aged.dll.L*(1+age*9);            % L    ..  10L
    aged.dll.eEps = aged.dll.eEps*(1-age/2);      % eps  ..  0.5eps

    % Particle cracking at positive electrode.
    aged.pos.Rs = aged.pos.Rs*(1-age*0.9);        % Rs   ..  0.1Rs   
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

    models(k) = setCellParam(model0,aged);
end

% Save model.
save(fullfile('agemodel','cellLMO_AgeSeries.mat'),'models','ageVect');