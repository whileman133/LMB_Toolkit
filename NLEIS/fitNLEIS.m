% fitNLEIS.m
%
% Validate estimation of reaction symmetry factors by regressing the 
% NLEIS model to COMSOL simulation data.
%
% -- Changelog --
% 2023.05.30 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath('..');
TB.addpaths();

% Load simulated lab spectra.
simName = 'cellLMO-Lumped-MSMR-30mA-socSeries';
load(fullfile('simdata',[simName '.mat']));
Z2lab = zeros(length(simData.freq),length(simData.socPct));
for k = 1:length(simData.socSeries)
    data = processEIS(simData.socSeries(k), 'NumHarmonics',2);
    Z2lab(:,k) = [data.h2.Zcell];
end

% True cell model for use in verification.
cellModel = evalSetpoint(simData.cellModel,[],0.5,simData.TdegC+273.15);
cellModelFlat = rmfield(cellModel,'function');

% Build optimization model.
param.neg.alpha = fastopt.param;
param.pos.alpha = fastopt.param('len',length(cellModel.pos.alpha));
modelspec = fastopt.modelspec(param);

% Build lower bounds.
lb.neg.alpha = 0.1;
lb.pos.alpha = 0.1*ones(size(cellModel.pos.alpha));
lb = fastopt.pack(lb,modelspec);

% Build upper bounds.
ub.neg.alpha = 0.9;
ub.pos.alpha = 0.9*ones(size(cellModel.pos.alpha));
ub = fastopt.pack(ub,modelspec);

fitModel = fastopt.particleswarm( ...
    @(params)cost(params,simData.freq,Z2lab,simData.socPct,simData.TdegC,cellModelFlat), ...
    modelspec,lb,ub, ...
    'particleCount',100,'swarmIterations',100, ...
    'fminconFunEvals',1000);
fitData.fitModel = fitModel;
fitData.trueModel = cellModel;
fitData.Z2lab = Z2lab;
fitData.simData = simData;
fitData.modelspec = modelspec;
fitData.lb = lb;
fitData.ub = ub;

% Save results to disk.
fileName = fullfile( ...
    'fitdata', ...
    [simName '.mat'] ...
);
save(fileName,"fitData");


function J = cost(params,freq,Z2lab,socPct,TdegC,cellModel)
    persistent bestJ;
    if isempty(bestJ)
        bestJ = Inf;
    end

    % Update model parameters
    mod = updateModel(cellModel,params);

    % Compute model-predicted second-harmonic impedance.
    Z2sim = zeros(size(Z2lab));
    for k = 1:length(socPct)
        % Compute second-harmonic spectra.
        tf = tfLMB(1j*2*pi*freq,mod,'Calc22',true, ...
            'TdegC',TdegC,'socPct',socPct(k));
        Z2sim(:,k) = tf.h22.tfVcell().';    
    end

    % Compute cost as sum of squared normalized residuals.
    J = sum(((Z2lab-Z2sim)./Z2lab).^2,'all');

    if J<bestJ
        bestJ = J;
        figure(1); clf;
        plot(real(Z2sim),-imag(Z2sim),'b-'); hold on;
        plot(real(Z2lab),-imag(Z2lab),'r.');
        drawnow;
        fprintf('J=%.3e, pos.alpha=%.2f %.2f, neg.alpha=%.2f\n', ...
            J, params.pos.alpha(1), params.pos.alpha(2), params.neg.alpha);
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