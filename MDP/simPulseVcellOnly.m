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
I = 1;                      % Amplitude of Iapp pulse [A].
tmax = 6;                   % Time when simulation ends [s].
tp = 1;
Ts = 0.1;
depletion = [0.1 0.5 1 2 10];
suffix = 'ce0';
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
kn0 = cellModel.function.neg.k0([],TdegC+273.15);
alphan0 = cellModel.function.neg.alpha([],TdegC+273.15);
kp0 = cellModel.function.pos.k0([],TdegC+273.15);
alphap0 = cellModel.function.pos.alpha([],TdegC+273.15);
qep = qep0*depletion;
qes = qep0*depletion;
qed = qed0*depletion;
psi = psi0*depletion;
kn = kn0*depletion.^(1-alphan0);
kp = kp0.*depletion.^(1-alphap0);

tt = 0:Ts:tmax;
iapp = zeros(size(tt));
iapp(tt>=tp) = I;
simspec.time = tt;
simspec.Iapp = iapp;
simspec.SOC0 = socPct;
simspec.T = TdegC;
simspec.TSHIFT = 0;

data = [];
for k = length(qep):-1:1
    kpString = sprintf('%g;',kp(:,k));
    mod = cellModel;
    mod.function.pos.qe = eval(sprintf('@(x,T)(%g)',qep(k)));
    mod.function.sep.qe = eval(sprintf('@(x,T)(%g)',qes(k)));
    mod.function.DL.qe = eval(sprintf('@(x,T)(%g)',qed(k)));
    mod.function.const.psi = eval(sprintf('@(x,T)(%g)',psi(k)));
    mod.function.neg.k0 = eval(sprintf('@(x,T)(%g)',kn(k)));
    mod.function.pos.k0 = eval(sprintf('@(x,T)([%s])',kpString));
    modelCOMSOL = genFOM(mod);
    [m1,sim] = simFOM(modelCOMSOL,simspec,OptSimFOM);
    m1.save(sprintf('C:\\Users\\whileman\\Desktop\\TOOLBOX_LMB\\PULSE\\PulseModel%d.mph',k));
    data(k).time = tt;
    data(k).Iapp = iapp;
    data(k).Vcell = sim.Vcell;
end

simData.pulse = data;
simData.depletion = depletion;
simData.cellModel = cellModel;
simData.TdegC = TdegC;
simData.socPct = socPct;

% Save results to disk.
[~,modelName,~] = fileparts(cellFile);
fileName = fullfile( ...
    'simdata', ...
    sprintf('%s-%dpct-%dmA', ...
        modelName,round(socPct),round(I*1000)) ...
);
if isfield(OptSimFOM,'FixExchangeCurrent') && OptSimFOM.FixExchangeCurrent
    fileName = [fileName '-fixI0'];
end
if isfield(OptSimFOM,'VcellOnly') && OptSimFOM.VcellOnly
    fileName = [fileName '-VcellOnly'];
end
if ~isempty(suffix)
    fileName = [fileName '-' suffix];
end
save(fileName,"simData");