% plotEdgeResistance.m
%
% Compute the edge resistances corresponding to various SOC and
% pulse-current setpoints for a simulation cell
%
% -- Changelog --
% 2024.01.06 | Created | Wesley Hileman <whileman@uccs.edu>

clear; clc; close all;
addpath(fullfile('..','..'));
TB.addpaths;

% Constants
cellModelFile = 'cellLMO-P2DM';
socPct = linspace(0, 100, 20);
socRctPct = linspace(0, 100, 100);
iappCrate = linspace(-5, 5, 20);
TdegC = 25;
plotdir = fullfile('plots',cellModelFile);

% Load cell model.
model = loadCellModel(cellModelFile);
model = convertCellModel(model,'WRM');
p.neg.Rdl = 0.4;
p.pos.Rdl = 0.1;
model = setCellParam(model,p);
[Q,theta0,theta100] = getCellParams( ...
    model,'const.Q,pos.theta0,pos.theta100','Output','list');
iapp = iappCrate*Q;

% Compute edge resistance.
Rinf = getPulseResistance(model,socPct,iapp,TdegC);

% Compute porous-electrode Rct for reference.
llpm = convertCellModel(model,'LLPM');
electrode = MSMR(llpm.function.pos);
theta = theta0 + (socRctPct/100)*(theta100-theta0);
RctData = electrode.Rct(llpm.function.pos,'theta',theta,'TdegC',TdegC);

% Plotting ----------------------------------------------------------------
if ~isfolder(plotdir)
    mkdir(plotdir);
end

figure; 
surf(socPct,iappCrate,Rinf); 
xlabel('State of Charge [%]'); 
ylabel('Magnitude [C rate]');
zlabel('R_\infty [\Omega]');
title('Edge Resistance vs. SOC and Current Magnitude')
thesisFormat3d;
print(fullfile(plotdir,'Rinf'),'-depsc');
print(fullfile(plotdir,'Rinf'),'-dpng');

figure;
plot(socRctPct,RctData.i0); hold on;
plot(socRctPct,RctData.i0j,'--');
xlabel('State of Charge, z [%]');
ylabel('Exchange Current, i_0^p [A]');
title('MSMR Exchange Current vs. SOC');
legend('\Sigma','j=1','j=2','Location','best');
thesisFormat;
print(fullfile(plotdir,'i0'),'-depsc');
print(fullfile(plotdir,'i0'),'-dpng');