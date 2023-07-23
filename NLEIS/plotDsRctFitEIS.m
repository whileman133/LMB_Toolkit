% plotDsRctFitEIS.m
%
% Plot regressed EIS model against lab data.

clear; close all; clc;
addpath(fullfile('..','..'));
TB.addpaths;

files = {
    fullfile('labfitdata','EIS-Cell395534_2-16degC.mat')
    fullfile('labfitdata','EIS-Cell395534-26degC.mat')
    fullfile('labfitdata','EIS-Cell395524-42degC.mat')
};
plotdir = fullfile('plots',sprintf('LAB-AllCells'));
if ~isfolder(plotdir)
    mkdir(plotdir);
end

soc = linspace(0,1,100);
socPct = soc*100;
for k = length(files):-1:1
    fitData = load(files{k});
    TdegC = fitData.TdegC;
    theta0 = fitData.values.pos.theta0;
    theta100 = fitData.values.pos.theta100;
    theta = theta0 + soc*(theta100-theta0);
    electrode = MSMR(fitData.values.pos);
    values(k) = fastopt.unpack( ...
        fastopt.pack(fitData.values,fitData.modelspec),fitData.modelspec, ...
        'flat',true,'sparse',true);
    lb(k) = fastopt.unpack( ...
        fastopt.pack(fitData.lb,fitData.modelspec),fitData.modelspec, ...
        'flat',true,'sparse',true);
    ub(k) = fastopt.unpack( ...
        fastopt.pack(fitData.ub,fitData.modelspec),fitData.modelspec, ...
        'flat',true,'sparse',true);
    ctData(k) = electrode.Rct(fitData.values.pos,'TdegC',TdegC,'theta',theta);
    dsData(k) = electrode.Ds(fitData.values.pos,'TdegC',TdegC,'theta',theta);
    socSpline(:,k) = 100*(fitData.values.pos.k0SplineTheta-theta0)/(theta100-theta0);
    RctSpline(:,k) = 1./fitData.values.pos.k0Spline./ctData(k).f;
    DsSpline(:,k) = fitData.values.pos.DsSpline;
    cellName = fitData.arg.labSpectra.cellName;
    cellName = strsplit(cellName,'_');
    cellName = cellName{1};
    labels{k} = sprintf('T=%.0f \\circC (#%s)',TdegC,cellName);
end

% Plot Rct.
% figure; colororder(cool(length(files)));
% semilogy(socPct,[ctData.Rct]); hold on;
% semilogy(socSpline,RctSpline,'o');
% ylabel('Charge-transfer resistance, R_{ct} [\Omega]');
% xlabel('Cell state-of-charge [%]');
% legend(labels{:},'Location','best');
% title(sprintf('Log-Spline Charge-Transfer'));
% thesisFormat;
% print(fullfile(plotdir,'Rct-soc'),'-depsc');
% print(fullfile(plotdir,'Rct-soc'),'-dpng');
% 
% % Plot Ds.
% figure; colororder(cool(length(files)));
% semilogy(socPct,[dsData.Ds]); hold on;
% semilogy(socSpline,DsSpline,'o');
% ylabel('Solid diffusion coefficient, D_s [s^{-1}]');
% xlabel('Cell state-of-charge [%]');
% legend(labels{:},'Location','best');
% title(sprintf('Log-Spline Solid Diffusion'));
% thesisFormat;
% print(fullfile(plotdir,'Ds-soc'),'-depsc');
% print(fullfile(plotdir,'Ds-soc'),'-dpng');

tabValues = struct2table(values);
tabValues.Label = labels(:);
tabValues = [tabValues(:,end) tabValues(:,1:end-1)];
tabLB = struct2table(lb);
tabLB.Label = labels(:);
tabLB = [tabLB(:,end) tabLB(:,1:end-1)];
tabUB = struct2table(ub);
tabUB.Label = labels(:);
tabUB = [tabUB(:,end) tabUB(:,1:end-1)];

writetable(tabValues,fullfile('labfitdata-xlsx','LinearEIS-AllCells.xlsx'),"Sheet",'values');
writetable(tabLB,fullfile('labfitdata-xlsx','LinearEIS-AllCells.xlsx'),"Sheet",'lb');
writetable(tabUB,fullfile('labfitdata-xlsx','LinearEIS-AllCells.xlsx'),"Sheet",'ub');