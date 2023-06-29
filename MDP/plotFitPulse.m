% plotFitPulse.m

clear; close all; clc;
if ~exist('TOOLBOX_LMB','dir')
    % bootstrap the toolbox
    addpath('..');
    TB.addPaths();
end

fileName = 'cellLMO-Lumped-MSMR-50pct-100mA-25degC-fixPsi';
load(fullfile('fitdata',[fileName '.mat']));
plotdir = fullfile('plots',fileName);
if ~isfolder(plotdir)
    mkdir(plotdir);
end

% Print estimation error.
comp = compareModels(fitData.trueModel,fitData.fitModel);
writetable(comp,fullfile('xlsx',[fileName '.xlsx']));

% Plot iapp(t).
figure;
plot(fitData.time,fitData.iapp/fitData.trueModel.const.Q);
xlabel('Time, t [s]');
ylabel('i_{app}(t) [C-rate]');
title('Charge-Neutral Pulse: i_{app}(t)');
thesisFormat([0.3 0 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,'iapp.eps'));
exportgraphics(gcf,fullfile(plotdir,'iapp.png'));

% Plot vcell(t) in both intervals;
for k = 1:length(fitData.regionsOfInterest)
    roi = fitData.regionsOfInterest{k};
    figure;
    plot(fitData.time(roi),fitData.vcell(roi));
    xlabel('Time, t [s]');
    ylabel('v_{cell}(t) [V]');
    title(sprintf('Charge-Neutral Pulse: v_{cell}(t) (Part %d)',k));
    thesisFormat([0.3 0 0.1 0.1]);
    exportgraphics(gcf,fullfile(plotdir,['vcell' num2str(k) '.eps']));
    exportgraphics(gcf,fullfile(plotdir,['vcell' num2str(k) '.png']));
end


function comparisonTable = compareModels(trueModel,fitModel)
    %COMPAREMODELS Compare true model parameter values to estimates.
    % Recursive implementation.

    cursor = 1;
    ParameterName = {};
    Truth = [];
    Estimate = [];
    PercentError = [];
    compare(trueModel,fitModel);
    ParameterName = ParameterName(:);
    Truth = Truth(:);
    Estimate = Estimate(:);
    PercentError = PercentError(:);
    comparisonTable = table(ParameterName,Truth,Estimate,PercentError);

    function compare(truth,est,prefix)
        if exist('prefix','var')
            prefix = [prefix '.'];
        else
            prefix = ''; 
        end
        paramNames = fieldnames(est);
        for k = 1:length(paramNames)
            paramName = paramNames{k};
            fullName = [prefix paramName];
            estValue = est.(paramName);
            trueValue = truth.(paramName);
            if isstruct(estValue)
                compare(trueValue,estValue,fullName);
            else
                pctError = 100*(trueValue-estValue)/trueValue;
                ParameterName{cursor} = fullName;
                Truth(cursor) = trueValue;
                Estimate(cursor) = estValue;
                PercentError(cursor) = pctError;
                cursor = cursor + 1;
            end % else
        end % for
    end % compare()
end % compareModels()