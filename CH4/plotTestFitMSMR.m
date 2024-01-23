% plotTestFitMSMR.m

clear; close all; clc;
addpath('..');
TB.addpaths;
load('testFitMSMR.mat');

truth = truth.toStruct();
clear stats;
for k = length(data):-1:1
    [err,rmse] = getPercentError(data(k).est,truth);
    stats(k).pcterr = err;
    stats(k).rmse = rmse;
    stats(k).Jcost = data(k).Jcost; 
end

% Vector of RMSE for each estimate (RMSE in percent).
rmse = [stats.rmse];

% Calculate success rate.
success = rmse<=0.1;
succrate = 100*sum(success)/length(rmse);

% Fetch best estimate.
[~,indbestest] = min(rmse);
bestest = stats(indbestest);


function [err, rmse] = getPercentError(estimate,truth)
%GETRMSE Calculate RMSE between estimates and true parameter values.

paramnames = fieldnames(truth);
toterr = 0;
totvar = 0;
for k = 1:length(paramnames)
    pname = paramnames{k};
    tru = truth.(pname);
    est = estimate.(pname);
    err.(pname) = 100*(tru-est)./tru;
    toterr = toterr + sum(err.(pname).^2);
    totvar = totvar + length(err.(pname));
end
rmse = sqrt(toterr/totvar);

end