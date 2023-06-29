% plotOCP.m
%
% Plot regressed MSMR models versus laboratory measurements.
%
% 2022.06.30 | Created |
% Wesley Hileman <whileman@uccs.edu>
% University of Colorado Colorado Springs

clear; close all; clc;

DATADIR = fullfile('..','modelfits');
PLOTDIR = fullfile('..','plots','modelfits');

loadOCP('fits','SionFresh_0C01-chg-J7','modelfits');
modelfits = TemperatureIndexedObjectWrapper(modelfits.name,modelfits.objects);
spreadsheetFile = fullfile(DATADIR,sprintf('fitMSMR-%s.xlsx',modelfits.name));

[temperatures, idxSort] = sort(modelfits.temperatures);
models = [modelfits.objects.model];
models = models(idxSort);

Uj0 = [models.Uj0];
Xj = [models.Xj];
Wj = [models.Wj];
J = models(1).J;
labelsGalleries = arrayfun(@(j)sprintf('j=%d',j),1:J,'UniformOutput',false);
labelsTemperatures = arrayfun(@(T)sprintf('%.0f\\circC',T),temperatures,'UniformOutput',false);

% figure; colororder(parula(J));
% plot(temperatures,Uj0,'d:','MarkerFaceColor','auto');
% ylim([3.5 4.7]);
% xlabel('Temperature [\circC]'); ylabel('U_j^0 [V]'); legend(labelsGalleries);
% title('LMB Cell, U_j^0');
% thesisFormat([0.2 0.1 0.1 0.2]);
% todisk(gcf,fullfile(PLOTDIR,sprintf('fitMSMR-%s-Uj0',modelfits.name)));
% 
% figure; colororder(parula(J));
% plot(temperatures,Xj,'d:','MarkerFaceColor','auto');
% xlabel('Temperature [\circC]'); ylabel('X_j [-]'); legend(labelsGalleries);
% title('LMB Cell, X_j');
% thesisFormat([0.2 0.1 0.1 0.2]);
% todisk(gcf,fullfile(PLOTDIR,sprintf('fitMSMR-%s-Xj',modelfits.name)));
% 
% figure; colororder(parula(models(1).J));
% plot(temperatures,Wj,'d:','MarkerFaceColor','auto');
% xlabel('Temperature [\circC]'); ylabel('\omega_j [-]'); legend(labelsGalleries);
% title('LMB Cell, \omega_j');
% thesisFormat([0.2 0.1 0.1 0.2]);
% todisk(gcf,fullfile(PLOTDIR,sprintf('fitMSMR-%s-Wj',modelfits.name)));
% 
% % Compare models at fit temperatures.
% figure; colororder(cool(length(temperatures)));
% for msmrfit = modelfits.objects
%     [Uocp, dzdv, Z, xj] = msmrfit.model.getOCP('volt',3.2,5,100000,msmrfit.temp);
%     subplot(121); plot(Z,Uocp); hold on;
%     subplot(122); plot(Uocp,abs(dzdv)); hold on;
% end
% subplot(121); 
% xlim([0 1]);
% xlabel('Absolute Stiochiometry, \theta');
% ylabel('Potential versus Li/Li+ [V]');
% legend(labelsTemperatures);
% title('Regressed OCP Models');
% subplot(122); 
% xlabel('Potential versus Li/Li+ [V]');
% ylabel('Differential Capacity, |d\theta/dU| [V^{-1}]');
% title('Differential Capacity');
% xlim([3.2 5]); ylim([0 3]);
% thesisFormat([0.2 0.1 0.1 0.1]);
% todisk(gcf,fullfile(PLOTDIR,sprintf('fitMSMR-%s-compare-models',modelfits.name)));
% 
% % Compare models at T=25degC (should be very similar if parameters do not
% % vary with temperature).
% figure; colororder(cool(length(temperatures)));
% for msmrfit = modelfits.objects
%     [Uocp, dzdv, Z, xj] = msmrfit.model.getOCP('volt',3.2,5,100000,25);
%     subplot(121); plot(Z,Uocp); hold on;
%     subplot(122); plot(Uocp,abs(dzdv)); hold on;
% end
% subplot(121); 
% xlim([0 1]);
% xlabel('Absolute Stiochiometry, \theta');
% ylabel('Potential versus Li/Li+ [V]');
% legend(labelsTemperatures);
% title('Regressed OCP Models (Evaluated at T=25\circC)');
% subplot(122); 
% xlabel('Potential versus Li/Li+ [V]');
% ylabel('Differential Capacity, |d\theta/dU| [V^{-1}]');
% title('Differential Capacity');
% xlim([3.2 5]); ylim([0 3]);
% thesisFormat([0.2 0.1 0.1 0.1]);
% todisk(gcf,fullfile(PLOTDIR,sprintf('fitMSMR-%s-compare-models-25C',modelfits.name)));

% for msmrfit = modelfits.objects
%     % Plot fit.
%     MSMR.compare(3,5,2,msmrfit.temp, ...
%         sprintf('Lab+XPD (%s %.0f\\circC)',msmrfit.name,msmrfit.temp),msmrfit.ocpest, ...
%         'Model Fit',msmrfit.model);
%     thesisFormat([0.2 0.1 0.1 0.1]);
%     todisk(gcf,fullfile(PLOTDIR,sprintf('fitMSMR-%s-%.0fC',modelfits.name,msmrfit.temp)));
% 
%     % Write model data to spreadsheet.
%     data = [
%         msmrfit.model.Uj0(:)';
%         msmrfit.model.Xj(:)';
%         msmrfit.model.Wj(:)';
%         msmrfit.model.zmin, msmrfit.model.zmax, nan(1,msmrfit.model.J-2);
%     ];
%     writematrix(data, spreadsheetFile,'Sheet',sprintf('T=%.0fdegC',msmrfit.temp));
% end

% Plot diagonal estimates (augmented with XPD data, so we need to
% normalize the relative stoichiometry of the laboratory data for 
% consistency).
figure;
colors = cool(length(modelfits.objects));
load(fullfile('..','processed','UocpXPD.mat'));
minZ_XPD = min(x);
for k = 1:length(modelfits.objects(idxSort))
    msmrfit = modelfits.objects(k);
    estimate = msmrfit.ocpest;
    test = estimate.ocptest;
    disZ = 0.18+test.disZ*(1-0.18);
    %disZ = (disZ-minZ_XPD)/(1-minZ_XPD);
    chgZ = 0.18+test.chgZ*(1-0.18);
    %chgZ = (chgZ-minZ_XPD)/(1-minZ_XPD);
    plot(minZ_XPD+estimate.Z*(1-minZ_XPD),estimate.V,'Color',colors(k,:),'DisplayName',labelsTemperatures{k}); hold on;
    plot(disZ,test.disV,':','Color',colors(k,:),'HandleVisibility','off','LineWidth',3);
    plot(chgZ,test.chgV,':','Color',colors(k,:),'HandleVisibility','off','LineWidth',3);
end
yline(4.3,'HandleVisibility','off'); xline(0.18,'HandleVisibility','off');
ylim([3.2 4.5]); xlim([0 1]);
xlabel('Absolute Lithiation of NMC Electrode, \theta');
ylabel('Potential versus Li/Li+ [V]');
title('Diagonal-Averaged OCP Estimates');
legend show;
thesisFormat([0.2 0.1 0.1 0.1]);

% Inset plot.
axes('position',[0.19 0.22 0.35 0.29]); box on;
for k = 1:length(modelfits.objects(idxSort))
    msmrfit = modelfits.objects(k);
    estimate = msmrfit.ocpest;
    test = estimate.ocptest;
    disZ = 0.18+test.disZ*(1-0.18);
    %disZ = (disZ-minZ_XPD)/(1-minZ_XPD);
    chgZ = 0.18+test.chgZ*(1-0.18);
    %chgZ = (chgZ-minZ_XPD)/(1-minZ_XPD);
    plot(minZ_XPD+estimate.Z*(1-minZ_XPD),estimate.V,'Color',colors(k,:),'DisplayName',labelsTemperatures{k}); hold on;
    plot(disZ,test.disV,':','Color',colors(k,:),'HandleVisibility','off','LineWidth',3);
    plot(chgZ,test.chgV,':','Color',colors(k,:),'HandleVisibility','off','LineWidth',3);
end
xlim(minZ_XPD+[0.08 0.12]);
ylim([4.185 4.275]);
todisk(gcf,fullfile(PLOTDIR,'OCPestimates'));