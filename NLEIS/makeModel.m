% makeModel.m
%
% Build state-space model for the system to which we'll apply MPC.
%
% -- Changelog --
% 2023.05.02 | Created | Wesley Hileman

clear; close all; clc;

% Generate model.
A = [1.1 2.0; 0 0.95];
B = [0; 0.0787];
C = [-1 1];
sys = ss(A,B,C,0);
poleOL = eig(A);
zeroOL = zero(sys);
save('MODEL.mat','sys','poleOL','zeroOL');

% Plot open-loop poles/zeros!
figure;
plot(complex(poleOL),'rx'); hold on;
plot(complex(zeroOL),'bo');
plot(complex(0+0j),'k+'); % origin
xline(0);
yline(0);
title('Open-Loop Pole/Zero Map');
xlabel('Re(s)');
ylabel('Im(s)');
setAxesNyquist;
thesisFormat;

% Plot open-loop step response.
