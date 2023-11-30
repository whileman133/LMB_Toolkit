% makeAgedModelLMO.m
%
% Make lumped-parameter models describing a full LMB cell over the course
% of cycle aging. LMO porous electrode.
%
% -- Changelog --
% 2023.09.28 | Eliminate reaction symmetry from age series | Wes H.
% 2023.09.14 | Sepatate LAM and LLI processes | Wes H.
% 2023.09.07 | Update for gen2 toolkit | Wes H.
% 2023.06.06 | Created | Wesley Hileman <whileman@uccs.edu>
%
% -- References --
% We use parameter values from this literature:
% [1] Shanshan Xu et al 2019 J. Electrochem. Soc. 166 A3456
% [2] Daniel R. Baker and Mark W. Verbrugge 2021 J. Electrochem. Soc. 168 050526
% [3] Dongliang Lu et al 2022 J. Electrochem. Soc. 169 080504

clear; close all; clc;
addpath('..');
TB.addpaths;

% Normalized variables that describe state-of-age of the cell.
lamVect = linspace(0,1,3);    % degree of loss of active material (LAM)
lliVect = linspace(0,1,3);    % degree of loss of lithium inventory (LLI)

% Model of fresh cell.
model0 = loadCellModel('cellLMO-P2DM');

clear ageArray;
for i = length(lamVect):-1:1
    for j = length(lliVect):-1:1
        % Fetch degrees of LAM and LLI.
        lam = lamVect(i);
        lli = lliVect(j);

        % Initialize struct of aged parameters.
        aged = getCellParams(model0);
    
        % 1. Loss of Active Material (LAM) at porous electrode --------
        % Particle cracking, decrease in diffusivity
        aged.pos.Rs = aged.pos.Rs*(1-lam*0.9);        % Rs   ..  0.1Rs  
        aged.pos.sEps = aged.pos.sEps*(1-lam*0.3);    % eps  ..  0.7eps
        aged.pos.Dsref = aged.pos.Dsref*(1-lam*0.8);  % Ds   ..  0.2Ds
        % Phase change in porous-electrode lattice.
        % NOTE: need to update theta0 and theta100 to obtain same 
        % voltage limits Vmin and Vmax!
        aged.pos.U0 = [
            aged.pos.U0(1)*(1+lam*0.02)   % 0%  ..  +2% drift
            aged.pos.U0(2)*(1+lam*0.01)   % 0%  ..  +1% drift
        ];
        aged.pos.omega = aged.pos.omega*(1+lam*0.1);  % 0%  ..  +10% drift
        ocpData = MSMR(aged.pos).ocp('voltage',[aged.const.Vmin aged.const.Vmax]);
        aged.pos.theta0 = ocpData.theta(1);   % theta @ Uocp(new)=Vmin
        aged.pos.theta100 = ocpData.theta(2); % theta @ Uocp(new)=Vmax
    
        % 2. Loss of Lithium Inventory (LLI) --------------------------
        % Electrolyte depletion.
        aged.const.ce0 = aged.const.ce0*(1-lli*0.5);    % ce0  ..  0.5ce0
        aged.const.kappa = aged.const.kappa*(1-lli*0.7);% k    ..  0.3k
        aged.const.De = aged.const.De*(1-lli*0.7);      % De   ..  0.3De
        % Dendrites and dead-lithium accumulation at negative electrode.
        aged.neg.gamma = aged.neg.gamma*(1+lli*2);      % gam  ..  3gam
        aged.neg.Rf = aged.neg.Rf*(1+lli*4);            % Rf   ..  5Rf
        aged.dll.L = aged.dll.L*(1+lli*9);              % L    ..  10L
        aged.dll.eEps = aged.dll.eEps*(1-lli/2);        % eps  ..  0.5eps
        % SEI film formation at positive electrode.
        aged.pos.Rf = aged.pos.Rf*(1+lli*4);            % Rf   ..  5Rf
    
        % Store aged model parameters.
        ageArray(i,j).model = setCellParam(model0,aged);
        ageArray(i,j).lam = lam;
        ageArray(i,j).lli = lli;
    end % for
end % for

% Save model.
save( ...
    fullfile('agemodel','cellLMO_AgeArray.mat'), ...
    'ageArray','lamVect','lliVect');