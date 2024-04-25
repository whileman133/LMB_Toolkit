% plotCompareROM_FOM.m

clear; close all; clc;
addpath('..');
TB.addpaths;

% Constants ---------------------------------------------------------------
simFile = 'cellLMO-P2DM_defaultHRA';
% Variables and x-locations to plot vs time.
vars.PhieTilde = [1 2 3];
vars.Phise = [2 3];
vars.Thetae = [0 1 2 3];
vars.Thetass = [2 3];
vars.Ifdl = [2 2.5 3];
vars.If = [2 2.5 3];
%vars.Eta = [1 2 3];
vars.Phis = [2 3];

% Collect simulation data -------------------------------------------------
simData = load(fullfile('SIM_FILES',[simFile '.mat']));
ROMout = simData.ROMout;
FOMout = simData.FOMout;
% Fetch time and cell voltage vectors.
timeHr = FOMout.time/60/60;
vcellROM = ROMout.Vcell;
vcellFOM = FOMout.Vcell;
% Fetch other variables.
varData = struct;
varNames = fieldnames(vars);
for k = 1:length(varNames)
    varName = varNames{k};
    xlocDesired = vars.(varName);
    varFOM = FOMout.(varName);
    varROM = ROMout.(varName);
    xlocFOM = FOMout.xLocs.(varName);
    xlocROM = ROMout.xLocs.(varName);
    [residFOM,indFOM] = min(abs(xlocFOM(:)-xlocDesired));
    [residROM,indROM] = min(abs(xlocROM(:)-xlocDesired));
    if any(residFOM>0.02)
        error("Invalid x-locs for FOM variable: %s",varName);
    end
    if any(residROM>0.01)
        error("Invalid x-locs for ROM variable: %s",varName);
    end
    varData.(varName).xloc = xlocDesired;
    varData.(varName).FOM = varFOM(:,indFOM);
    varData.(varName).ROM = varROM(:,indROM);
end % for vars

% Plotting ----------------------------------------------------------------
plotdir = fullfile('PLOTS',simFile);
if ~isfolder(plotdir)
    mkdir(plotdir);
end

figure;
plot(timeHr,vcellROM); hold on;
plot(timeHr,vcellFOM,':');
title('vcell-v-time');
xlabel('time-hr');
ylabel('vcell');
legend('HRA','COMSOL','Location','best');
thesisFormat;
print(fullfile(plotdir,'vcell'),'-depsc');
print(fullfile(plotdir,'vcell'),'-dpng');

varNames = fieldnames(varData);
for k = 1:length(varNames)
    varName = varNames{k};
    var = varData.(varName);
    colorsROM = winter(length(var.xloc));
    colorsFOM = cool(length(var.xloc));
    figure;
    for j = 1:length(var.xloc)
        xloc = var.xloc(j);
        plot(timeHr,var.ROM(:,j), ...
            'Color',colorsROM(j,:), ...
            'DisplayName',sprintf('x=%.0f HRA',xloc)); 
        hold on;
    end % for
    for j = 1:length(var.xloc)
        xloc = var.xloc(j);
        plot(timeHr,var.FOM(:,j),':', ...
            'Color',colorsFOM(j,:), ...
            'DisplayName',sprintf('x=%.1f COMSOL',xloc)); 
    end % for
    xlabel('time-hr');
    ylabel('var');
    title(varName);
    legend('Location','best','NumColumns',2);
    thesisFormat;
    print(fullfile(plotdir,varName),'-depsc');
    print(fullfile(plotdir,varName),'-dpng');
end % for varData

% RMSE Computation --------------------------------------------------------
vcellRMSE = sqrt(mean((vcellROM-vcellFOM).^2));
fprintf('RMSE:%-13s = %10.3f mV\n','Vcell',vcellRMSE*1000);
varNames = fieldnames(varData);
for k = 1:length(varNames)
    varName = varNames{k};
    var = varData.(varName);
    for j = 1:length(var.xloc)
        xloc = var.xloc(j);
        rmse = sqrt(mean((var.ROM(:,j)-var.FOM(:,j)).^2));
        fprintf('RMSE:%-10s(%3.1f) = %10.3f m\n',varName,xloc,rmse*1000);
    end % for
end