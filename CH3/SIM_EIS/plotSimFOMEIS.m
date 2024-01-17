% plotSimFOMEIS.m
%
% Plot results of medium-signal EIS simulation for full-order LMB cell.
%
% -- Changelog --
% 2023.12.17 | Make code DRY with variable loop | Wes H.
% 2023.04.05 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath(fullfile('..','..'));
TB.addpaths;

% Constants.
simName = 'cellLMO-P2DM-1mA-10pct';
plotdir = fullfile('plots',simName);
socPctPlot = 10;   % SOC setpoint to plot
vars.Thetae = struct('xpos',[0 1 2 2.5 3]);
vars.Thetass = struct('xpos',[2 2.5 3]);
vars.PhieTilde = struct('xpos',[1 2 2.5 3]); % =0 at x=0, so exlcude
vars.PhisTilde = struct('xpos',[2 2.5]); % =0 at x=3, so exclude
vars.Phise = struct('xpos',[2 2.5 3]);
vars.Eta = struct('xpos',[2 2.5 3]);
vars.If = struct('xpos',[2 2.5 3]);
vars.Ifdl = struct('xpos',[2 2.5 3]);

% Fetch and process simulation data.
% Also evalulate linear TFs of same variables at same x-locs for 
% comparison to COMSOL simulation.
seriesData = load(fullfile('simdata',[simName '.mat']));
seriesData = seriesData.simData;
[~,indSOC] = min(abs(seriesData.socPct-socPctPlot));
simData = seriesData.socSeries(indSOC);
modelName = simData.arg.cellModel.metadata.cell.name;
spectra = processEIS(simData, ...
    'NumHarmonics',2,'EvalLinTF',true,'NumTFFreqPoints',200);

% Genrate directory in which to place plots.
if ~isfolder(plotdir)
    mkdir(plotdir);
end


% Plot linear TFs ---------------------------------------------------------

% Plot linear impedance (Nyquist).
figure;
plot(real([spectra.tf.Zcell]),-imag([spectra.tf.Zcell])); hold on;
plot(real([spectra.lin.Zcell]),-imag([spectra.lin.Zcell]),'d');
title(sprintf('Zcell'));
xlabel('Re');
ylabel('-Im');
legend('TF','FOM','Location','best');
setAxesNyquist;
thesisFormat('LineMarkerSize',5);
ax = addInset([3.36 3.59],[2.9 0.2]);
setAxesNyquist('axes',ax,'xdata',[3.36 3.59],'ydata',[0 0.07]);
print(fullfile(plotdir,'Zcell-Nyq.eps'),'-depsc');
print(fullfile(plotdir,'Zcell-Nyq.png'),'-dpng');

% Plot linear TFs for internal variables (Nyqiust).
varnames = fieldnames(vars);
for v = 1:length(varnames)
    varname = varnames{v};
    vardata = vars.(varname);
    xposPlot = vardata.xpos;
    xposSim = spectra.xlocs.(varname);
    [~,indxx] = min(abs(xposPlot(:)'-xposSim(:)));
    frTF = [spectra.tf.(varname)];
    frSIM = [spectra.lin.(varname)];
    labels1 = arrayfun(@(x)sprintf('TF'),xposPlot, ...
        'UniformOutput',false);
    labels2 = arrayfun(@(x)sprintf('FOM $$\\tilde{x}=%.1f$$',x),xposPlot, ...
        'UniformOutput',false);
    colors = cool(length(xposPlot));

    figure;
    for k = 1:length(xposPlot)
        plot(real(frTF(indxx(k),:)),-imag(frTF(indxx(k),:)), ...
            'Color',colors(k,:)); 
        hold on;
    end
    for k = 1:length(xposPlot)
        plot(real(frSIM(indxx(k),:)),-imag(frSIM(indxx(k),:)),'d', ...
            'Color',colors(k,:),'MarkerFaceColor','auto');
    end
    title(varname);
    xlabel('Re');
    ylabel('-Im');
    legend(labels1{:},labels2{:},'NumColumns',2,'Location','best', ...
        'Interpreter','latex');
    setAxesNyquist;
    thesisFormat('LineMarkerSize',4);
    if strcmpi(varname,'Phise')
        ax = addInset([-0.2 -0.04],[-0.75 -0.8],'YSpan',[-0.05 0]);
        setAxesNyquist('axes',ax,'xdata',[-0.2 -0.04],'ydata',[-0.05 0]);
    end
    if strcmpi(varname,'PhieTilde')
        ax = addInset([-0.013 -0.01],[-0.08 -0.04]);
        setAxesNyquist('axes',ax,'xdata',[-0.013 -0.01],'ydata',[-8e-4 0]);
    end
    print(fullfile(plotdir,[varname '-Nyq.eps']),'-depsc');
    print(fullfile(plotdir,[varname '-Nyq.png']),'-dpng');
end