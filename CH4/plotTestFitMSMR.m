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
Jcost = [stats.Jcost];

% Calculate success rate.
success = rmse<=0.1;
succrate = 100*sum(success)/length(rmse);

% Fetch best estimate.
[~,indbest] = min(rmse);
beststats = stats(indbest);
bestest = data(indbest);

% Fetch worst estimate.
[~,indworst] = max(rmse);
worststats = stats(indworst);
worstest = data(indworst);

% Plot best/worst models.
bm = MSMR(bestest.est);    bmOCP = bm.ocp('npoints',100);
wm = MSMR(worstest.est);   wmOCP = wm.ocp('npoints',100);
tm = MSMR(truth);          tmOCP = tm.ocp('npoints',100);

figure;
subplot(121);
plot(tmOCP.theta,tmOCP.Uocp,'k:'); hold on;
plot(bmOCP.theta,bmOCP.Uocp,'m--');
xlabel('theta');
ylabel('Uocp');
title('OCP');
legend('Truth','Model');
subplot(122);
plot(abs(1./tmOCP.dUocp),tmOCP.Uocp,'k:'); hold on;
plot(abs(1./bmOCP.dUocp),bmOCP.Uocp,'m--');
xlabel('absdZ');
ylabel('Uocp');
title('DC');
thesisFormat;
print(fullfile('plots','fitmsmrBEST'),'-depsc');
print(fullfile('plots','fitmsmrBEST'),'-dpng');

figure;
subplot(121);
plot(tmOCP.theta,tmOCP.Uocp,'k:'); hold on;
plot(wmOCP.theta,wmOCP.Uocp,'m--');
xlabel('theta');
ylabel('Uocp');
title('OCP');
legend('Truth','Model');
subplot(122);
plot(abs(1./tmOCP.dUocp),tmOCP.Uocp,'k:'); hold on;
plot(abs(1./wmOCP.dUocp),wmOCP.Uocp,'m--');
xlabel('absdZ');
ylabel('Uocp');
title('DC');
thesisFormat;
print(fullfile('plots','fitmsmrWORST'),'-depsc');
print(fullfile('plots','fitmsmrWORST'),'-dpng');

figure;
histogram(log10([stats.rmse]),25,'Normalization','percentage','FaceAlpha',1);
xline(-1,'r--');
xlabel('log10pcterror');
ylabel('relfreq');
title('distRMSE');
text(-1.02,17,'emax', ...
    'VerticalAlignment', 'top', 'HorizontalAlignment','right', ...
    'FontSize',16);
thesisFormat;
print(fullfile('plots','fitmsmrDIST'),'-depsc');
print(fullfile('plots','fitmsmrDIST'),'-dpng');

figure;
scatter(log10(rmse),log10(Jcost),100);
xlabel('log10pcterror');
ylabel('log10Jcost');
title('cost-v-pcterror');
thesisFormat;
print(fullfile('plots','fitmsmrCOSTVERR'),'-depsc');
print(fullfile('plots','fitmsmrCOSTVERR'),'-dpng');


% Print results.
printMSMR('Truth',truth);
printMSMR('Best Fit Model',bm.toStruct());
printMSMR('Percent Error',getPercentError(bm.toStruct(),truth));
printMSMR('Worst Fit Model',wm.toStruct());
printMSMR('Percent Error',getPercentError(wm.toStruct(),truth));


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

function printMSMR(name,params)
%PRINTMSMR Print MSMR parameters to the console.

fprintf("%s\n",upper(name));
fprintf("  U0      : ");
fprintf("%-10.3f",params.U0);
fprintf("\n");
fprintf("  X       : ");
fprintf("%-10.3f",params.X);
fprintf("\n");
fprintf("  omega   : ");
fprintf("%-10.3f",params.omega);
fprintf("\n");
fprintf("  thetamin: ");
fprintf("%-10.3f",params.thetamin);
fprintf("\n");
fprintf("  thetamax: ");
fprintf("%-10.3f",params.thetamax);
fprintf("\n");

end