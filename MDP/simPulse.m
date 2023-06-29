% simPulse.m
%
% Simulate medium-duration pulse response of the full-order LMB model 
% in COMSOL for several values of the salt inventory and save results 
% to disk.
%
% -- Changelog --
% 04.18.2023 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath(fullfile("..","UTILITY"));
addpath(fullfile("..","XLSX_CELLDEFS"));

% Constants.
cellFile = 'cellLMO-Lumped-MSMR-modk0.xlsx';  % Name of cell parameters spreadsheet.
socPct = 50;                % Cell SOC setpoint [%].
TdegC = 25;                 % Cell temperature [degC].
I = 0.1;                    % Amplitude of Iapp pulse [A].
tp1 = 1;                    % Time when iapp discharge pulse starts [s].
tp2 = 31;                   % Time when iapp charge pulse starts [s].
tmax = 61;                  % Time when simulation ends [s].
depletion = [0.1 0.5 1 2 10];
suffix = 'psikD';
% Structure of additional options to pass to simFOM.
OptSimFOM.FixExchangeCurrent = false;
OptSimFOM.VcellOnly = true;

% Load cell parameters and run EIS simulation in COMSOL.
cellModel = loadCellParams(cellFile);

% Porous-electrode salt-inventories to try.
qep0 = cellModel.function.pos.qe([],TdegC+273.15);
qes0 = cellModel.function.sep.qe([],TdegC+273.15);
qed0 = cellModel.function.DL.qe([],TdegC+273.15);
psi0 = cellModel.function.const.psi([],TdegC+273.15);
kappaD0 = cellModel.function.const.kD([],TdegC+273.15);
kn0 = cellModel.function.neg.k0([],TdegC+273.15);
alphan0 = cellModel.function.neg.alpha([],TdegC+273.15);
kp0 = cellModel.function.pos.k0([],TdegC+273.15);
alphap0 = cellModel.function.pos.alpha([],TdegC+273.15);
qep = qep0*depletion;
qes = qes0*depletion;
qed = qed0*depletion;
psi = psi0*depletion;
kappaD = kappaD0*depletion;
kn = kn0*depletion.^(1-alphan0);
kp = kp0.*depletion.^(1-alphap0);

tt = linspace(0,tmax,1000);
iapp = zeros(size(tt));
iapp(tp1<tt&tt<=tp2) = I;
iapp(tp2<tt) = -I;
simspec.time = tt;
simspec.Iapp = iapp;
simspec.SOC0 = socPct;
simspec.T = TdegC;
simspec.TSHIFT = 0;

data = [];
for k = length(depletion):-1:1
    kpString = sprintf('%g;',kp(:,k));
    mod = cellModel;
%     mod.function.pos.qe = eval(sprintf('@(x,T)(%g)',qep(k)));
%     mod.function.sep.qe = eval(sprintf('@(x,T)(%g)',qes(k)));
%     mod.function.DL.qe = eval(sprintf('@(x,T)(%g)',qed(k)));
     mod.function.const.psi = eval(sprintf('@(x,T)(%g)',psi(k)));
     mod.function.const.kD = eval(sprintf('@(x,T)(%g)',kappaD(k)));
%     mod.function.neg.k0 = eval(sprintf('@(x,T)(%g)',kn(k)));
%     mod.function.pos.k0 = eval(sprintf('@(x,T)([%s])',kpString));
    modelCOMSOL = genFOM(mod);
    [~,sim] = simFOM(modelCOMSOL,simspec,OptSimFOM);
    data(k).Vcell = sim.Vcell;
end

simData.pulse = data;
simData.depletion = depletion;
simData.cellModel = cellModel;
simData.TdegC = TdegC;
simData.socPct = socPct;
simData.time = simspec.time;
simData.iapp = simspec.Iapp;

% Save results to disk.
[~,modelName,~] = fileparts(cellFile);
fileName = fullfile( ...
    'simdata', ...
    sprintf('%s-%dpct-%dmA-%ddegC', ...
        modelName,round(socPct),round(I*1000),round(TdegC)) ...
);
if isfield(OptSimFOM,'FixExchangeCurrent') && OptSimFOM.FixExchangeCurrent
    fileName = [fileName '-fixI0'];
end
if ~isempty(suffix)
    fileName = [fileName '-' suffix];
end
save(fileName,"simData");