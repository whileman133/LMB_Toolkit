% plotPRSurface.m
%
% Plot a pulse-resistance surface (collection of R0 estimates over SOC
% and C-rate) based on data collected in the laboratory.
%
% -- Changelog --
% 07.23.2023 | Updated for gen2 toolkit | Wesley Hileman
% 08.17.2022 | Updated for new schema | Wesley Hileman
% 06.02.2022 | Created |
% Wesley Hileman <whileman@uccs.edu>
% University of Colorado Colorado Springs

clearvars -except study; clc; close all;
addpath('..');
TB.addpaths;

% Name of pulse-resistance estimates to load.
STUDYNAME = 'SionFresh-P0-2tau-100us-ss';

% Temperature at which to plot the estimate (selects dataset).
T = 0; % or 0, 15, 25, 40 degC

% Parameter to plot.
PARAM.dispname = 'R_{0}';
PARAM.dispunit = '\Omega';
PARAM.getFromECM = @(ecm)ecm.R0;

% Plot 3-sigma error bounds if true. 
SHOWBOUNDS = false;

% Saves plots to disk if true.
SAVEPLOTS = false;

% Directory in which to save plots.
PLOTDIR = fullfile(TB.const.SDPROOT,'labdata','plots','surfaces');

% Plot tick locations for C-rate.
CRATETICKS = -5:5;

% Whether or not to limit the range of the plotted parameter.
LIMITPARAM = true;

paramClean = erase(PARAM.dispname,{'{','}','_'});

pr.load(STUDYNAME,'surfaces','varname','estimates');
est = estimates.getTest(T);
currents = est.iapp(:);
socs = est.z(:)*100;
[R, Rbound] = est.getParam(PARAM.getFromECM);

% Outliers make the surface plot messy; attempt to remove them.
idxOutl = isoutlier(R,'mean',2);
rsurf = R;
rsurf(idxOutl) = nan;
limR = [min(rsurf,[],'all') max(rsurf,[],'all')];
%limR = [0.15  0.25];

idxI = find(currents < -0.04);
idxI = flip(idxI);
idxI = idxI(1:2:length(idxI));
idxZ = 1:length(socs);

wh.setFigure; figure; colororder(cool(length(idxI)));
plot(socs(idxZ),R(idxI,idxZ),':d'); hold on;
if SHOWBOUNDS
    colororder(cool(length(idxI)));
    plot( ...
        [socs(idxZ); nan; socs(idxZ)], ...
        [R(idxI,idxZ)-Rbound(idxI,idxZ), nan(length(idxI),1), R(idxI,idxZ)+Rbound(idxI,idxZ)], ...
        '-.');
end
if LIMITPARAM
    ylim(limR);
end
names = cell.empty(length(idxI),0);
for k = 1:length(idxI)
    names{k} = sprintf('%.2fC', currents(idxI(k))/est.QAh);
end
legend(names);
xlabel('SOC [%]');
ylabel(sprintf('%s [%s]',PARAM.dispname,PARAM.dispunit));
title(sprintf('%s: %s (%d\\circC)',PARAM.dispname,est.name,est.temp));
thesisFormat();
if SAVEPLOTS
    todisk(gcf,fullfile(PLOTDIR,sprintf('%s-v-Ichg-%s-%ddegC',paramClean,est.name,est.temp)));
end

idxI = find(currents > 0.04);
idxI = idxI(1:2:length(idxI));
idxZ = 1:length(socs);

wh.setFigure; figure; colororder(cool(length(idxI)));
plot(socs(idxZ),R(idxI,idxZ),':d'); hold on;
if SHOWBOUNDS
    colororder(cool(length(idxI)));
    plot( ...
        [socs(idxZ); nan; socs(idxZ)], ...
        [R(idxI,idxZ)-Rbound(idxI,idxZ), nan(length(idxI),1), R(idxI,idxZ)+Rbound(idxI,idxZ)], ...
        '-.');
end
if LIMITPARAM
    ylim(limR);
end
names = cell.empty(length(idxI),0);
for k = 1:length(idxI)
    names{k} = sprintf('%.2fC', currents(idxI(k))/est.QAh);
end
legend(names);
xlabel('SOC [%]');
ylabel(sprintf('%s [%s]',PARAM.dispname,PARAM.dispunit));
title(sprintf('%s: %s (%d\\circC)',PARAM.dispname,est.name,est.temp));
thesisFormat();
if SAVEPLOTS
    todisk(gcf,fullfile(PLOTDIR,sprintf('%s-v-Idis-%s-%ddegC',paramClean,est.name,est.temp)));
end

idxI = 1:length(currents);
idxZ = 1:3:length(socs);

figure; colororder(cool(length(idxZ)));
plot(currents(idxI)/est.QAh,R(idxI,idxZ),':d'); hold on;
if SHOWBOUNDS
    colororder(cool(length(idxZ)));
    plot( ...
         [currents(idxI)/est.QAh; nan; currents(idxI)/est.QAh], ...
         [R(idxI,idxZ)-Rbound(idxI,idxZ); nan(1,length(idxZ)); R(idxI,idxZ)+Rbound(idxI,idxZ)], ...
         '-.');
end
if LIMITPARAM
    ylim(limR);
end
xticks(CRATETICKS);
names = cell.empty(length(idxZ),0);
for k = 1:length(idxZ)
    names{k} = sprintf('%.0f%%', socs(idxZ(k)));
end
legend(names);
xlabel('C-rate');
ylabel(sprintf('%s [%s]',PARAM.dispname,PARAM.dispunit));
title(sprintf('%s: %s (%d\\circC)',PARAM.dispname,est.name,est.temp));
thesisFormat();
if SAVEPLOTS
    todisk(gcf,fullfile(PLOTDIR,sprintf('%s-v-SOC-%s-%ddegC',paramClean,est.name,est.temp)));
end

figure;
surf(socs(:),currents(:)/est.QAh,rsurf);
if LIMITPARAM
    zlim(limR);
end
xlabel('SOC [%]');
ylabel('C-rate');
zlabel(sprintf('%s [%s]',PARAM.dispname,PARAM.dispunit));
title(sprintf('%s: %s (%d\\circC)',PARAM.dispname,est.name,est.temp));
thesisFormat3d([0.2 0.1 0.3 0.1]);
if SAVEPLOTS
    todisk(gcf,fullfile(PLOTDIR,sprintf('%s-3d-%s-%ddegC',paramClean,est.name,est.temp)));
end