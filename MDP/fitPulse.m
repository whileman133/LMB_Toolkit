% fitPulse.m
%
% Validate estimation of salt inventory parameters (qep and qes) and kD
% using a medium-duration pulse by regressing the full-order COMSOL model 
% to the simulated pulse response.

clear; close all; clc;
addpath(fullfile("..","UTILITY"));
simname = 'cellLMO-Lumped-MSMR-modk0-50pct-100mA-25degC-psikD';
savename = 'cellLMO-Lumped-MSMR-50pct-100mA-25degC-fixPsi';
load(fullfile('simData',[simname '.mat']));
pulse = simData.pulse(find(simData.depletion==1,1,'first'));
time = simData.time;
iapp = simData.iapp;
vcell = pulse.Vcell;
cellModel = evalSetpoint( ...
    simData.cellModel,[], ...
    simData.socPct/100,simData.TdegC+273.15);

% Build optimization model.
params.pos.qe = fastopt.param;
params.sep.qe = fastopt.param;
params.const.kD = fastopt.param;
params.const.psi = fastopt.param('fix',cellModel.const.psi);
modelspec = fastopt.modelspec(params);

% Build lower bounds.
lb.pos.qe = cellModel.pos.qe/50;
lb.sep.qe = cellModel.sep.qe/50;
lb.const.kD = cellModel.const.kD*50;  % kD is negative!
lb.const.psi = cellModel.const.psi/50;
lb = fastopt.pack(lb,modelspec);

% Build upper bounds.
ub.pos.qe = cellModel.pos.qe*50;
ub.sep.qe = cellModel.sep.qe*50;
ub.const.kD = cellModel.const.kD/50;  % kD is negative!
ub.const.psi = cellModel.const.psi*50;
ub = fastopt.pack(ub,modelspec);

% Determine regression interval; 0.1sec after onset of pulses
% to avoid CPE effects of double-layers!
dis = iapp>0;  
chg = iapp<0;
tdis = time(find(iapp>0,1,'first'));  % time at onset of discharge pulse
tchg = time(find(iapp<0,1,'first'));  % time at onset of charge pulse
reg = false(size(time));  % logical indicies to regression interval
regions = {};
if ~isempty(tdis)
    dis = tdis+0.1<=time&time<=tdis+25;  % logical indicies to discharge interval
    reg = reg | dis;
    regions{end+1} = dis;
end
if ~isempty(tchg)
    chg = tchg+0.1<=time&time<=tchg+25;  % logical indicies to charge interval
    reg = reg | chg;
    regions{end+1} = chg;
end

model = fastopt.particleswarm( ...
    @(params)cost(params,cellModel,time,iapp,vcell,reg,regions), ...
    modelspec,lb,ub, ...
    'particleCount',100,'swarmIterations',50, ...
    'fminconFunEvals',1000);
fitData.fitModel = model;
fitData.trueModel = cellModel;
fitData.time = time;
fitData.iapp = iapp;
fitData.vcell = vcell;
fitData.regressionInterval = reg;
fitData.regionsOfInterest = regions;
fitData.socPct = simData.socPct;
fitData.TdegC = simData.TdegC;
fitData.modelspec = modelspec;
fitData.lb = lb;
fitData.ub = ub;

% Save results to disk.
fileName = fullfile( ...
    'fitdata', ...
    [savename '.mat'] ...
);
save(fileName,"fitData");

function J = cost(params,cellModel,timelab,iapplab,vcelllab,reg,regions)
    persistent bestJ;
    if isempty(bestJ)
        bestJ = Inf;
    end

    % Update model parameters.
    mod = cellModel;
    mod.function = updateModel(cellModel.function,params);

    % Run COMSOL simulation.
    simspec.time = timelab;
    simspec.Iapp = iapplab;
    simspec.SOC0 = cellModel.const.soc*100;
    simspec.T = cellModel.const.T-273.15;
    simspec.TSHIFT = 0;
    modelCOMSOL = genFOM(mod,'DebugFlag',false);
    [~,sim] = simFOM(modelCOMSOL,simspec,'VcellOnly',true,'DebugFlag',false);
    vcellsim = sim.Vcell;

    % Compute cost as sum of squared residuals.
    J = sum((vcellsim(reg)-vcelllab(reg)).^2);

    if J<bestJ
        bestJ = J;
        for k = 1:length(regions)
            figure(k); clf;
            plot(timelab(regions{k}),vcelllab(regions{k}),'b.'); hold on;
            plot(timelab(regions{k}),vcellsim(regions{k}),'r.');
            drawnow;
        end
        fprintf('J=%.5e pos.qe=%.5e sep.qe=%.5e psi=%.5e\n',J, ...
            params.pos.qe, params.sep.qe, params.const.psi);
    end
end

function newModel = updateModel(model,params)
    %UPDATEMODEL Update LMB model parameters.
    newModel = model;
    paramnames = fieldnames(newModel);
    for k = 1:length(paramnames)
        paramname = paramnames{k};
        if ~isfield(params,paramname)
            continue;
        end
        if isstruct(newModel.(paramname))
            newModel.(paramname) = updateModel( ...
                newModel.(paramname),params.(paramname));
        elseif isa(newModel.(paramname),'function_handle')
            newModel.(paramname) = eval( ...
                sprintf('@(x,T)(%g)',params.(paramname)));
        else
            newModel.(paramname) = params.(paramname);
        end
    end
end