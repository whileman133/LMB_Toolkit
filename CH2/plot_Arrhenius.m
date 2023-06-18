% plotArrhenius.m
%
% Wes Hileman
% ECE 5710
% 11.11.2021
%
% Plot the Arrhenius scale factor versus temperature and activation
% energy.

clear; close all; clc;

Tref = 25 + 273;                    % reference temperature [K]
TT = linspace(0, 50, 100) + 273;    % temperature vector [K]
EEa = [4000, 20000, 30000];         % activation energies [J/mol]
R = 8.3145;                         % ideal gas constant [J/mol K]

figure;
for Ea = EEa
    arr = exp(Ea*(1/Tref - 1./TT)/R);
    plot(TT - 273, arr, 'DisplayName', sprintf('E_a = %d J/mol', Ea)); 
    hold on;
end
legend('Location', 'northwest');
xlabel('Temperature [\circC]');
ylabel('Arrhenius Scale Factor');
title('Arrhenius Scale Factor vs Temperature and E_a');
thesisFormat([0.1 0.1 0.1 0.1]);
exportgraphics(gcf,'ArrheniusScale.eps');
exportgraphics(gcf,'ArrheniusScale.png');
