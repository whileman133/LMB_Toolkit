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
simName = 'cellLMO-P2DM-5mA-10pct';
plotdir = fullfile('plots',simName);
socPctPlot = 10;   % SOC setpoint to plot
vars.Thetae = struct('xpos',[0 1 2 2.5 3]);
vars.Thetass = struct('xpos',[2 2.5 3]);
vars.PhieTilde = struct('xpos',[1 2 2.5 3]); % =0 at x=0, so exlcude
%vars.PhisTilde = struct('xpos',[2 2.5]); % =0 at x=3, so exclude
vars.Phise = struct('xpos',[2 2.5 3]);
vars.Eta = struct('xpos',[2 2.5 3]);
%vars.If = struct('xpos',[2 2.5 3]);
vars.Ifdl = struct('xpos',[2 2.5 3]);

% Fetch and process simulation data.
% Also evalulate second-harmonic phasors of same variables at same x-locs 
% for comparison to COMSOL simulation.
seriesData = load(fullfile('simdata',[simName '.mat']));
seriesData = seriesData.simData;
[~,indSOC] = min(abs(seriesData.socPct-socPctPlot));
simData = seriesData.socSeries(indSOC);
cellModel = simData.arg.cellModel;
modelName = cellModel.metadata.cell.name;
spectra = processEIS(simData, ...
    'NumHarmonics',2,'EvalLinTF',false,'NumTFFreqPoints',200);
freqSim = spectra.freq;
freq = logspace(log10(min(freqSim)),log10(max(freqSim)),100);
tf = tfLMB(1j*2*pi*freq,cellModel,'Calc22',true, ...
    'TdegC',simData.arg.TdegC,'socPct',simData.arg.socPct);

% Genrate directory in which to place plots.
if ~isfolder(plotdir)
    mkdir(plotdir);
end


% Plot linear TFs ---------------------------------------------------------

% Plot linear impedance (Nyquist).
tfZcell = tf.h22.tfZcell().';
figure;
plot(real(tfZcell),-imag(tfZcell)); hold on;
plot(real([spectra.h2.Zcell]),-imag([spectra.h2.Zcell]),'d');
title('Zcell');
xlabel('Re');
ylabel('-Im');
legend('H2','FOM','Location','best');
setAxesNyquist;
thesisFormat('LineMarkerSize',5);
%ax = addInset([3.36 3.59],[2.9 0.2]);
%setAxesNyquist('axes',ax,'xdata',[3.36 3.59],'ydata',[0 0.07]);
print(fullfile(plotdir,'Zcell-Nyq.eps'),'-depsc');
print(fullfile(plotdir,'Zcell-Nyq.png'),'-dpng');

% Plot linear TFs for internal variables (Nyqiust and Bode).
varnames = fieldnames(vars);
for v = 1:length(varnames)
    varname = varnames{v};
    vardata = vars.(varname);
    xposPlot = vardata.xpos;
    xposSim = spectra.xlocs.(varname);
    [~,indxx] = min(abs(xposPlot(:)'-xposSim(:)));
    frTF = tf.h22.(['tf' varname])(xposSim);
    frTF = frTF.'; % !!! important to transpose
    frSIM = [spectra.h2.(varname)];
    labels1 = arrayfun(@(x)sprintf('H2'),xposPlot, ...
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
    print(fullfile(plotdir,[varname '-Nyq.eps']),'-depsc');
    print(fullfile(plotdir,[varname '-Nyq.png']),'-dpng');

%     figure;
%     for k = 1:length(xposPlot)
%         loglog(freq,abs(frTF(indxx(k),:)),'Color',colors(k,:)); 
%         hold on;
%     end
%     for k = 1:length(xposPlot)
%         loglog(freqSim,abs(frSIM(indxx(k),:)),'d','Color',colors(k,:)); 
%     end
%     title(varname);
%     xlabel('Freq');
%     ylabel('Mag');
%     legend(labels1{:},labels2{:},'NumColumns',2,'Location','best', ...
%         'Interpreter','latex');
%     thesisFormat('LineMarkerSize',4);
%     print(fullfile(plotdir,[varname '-BodeMag.eps']),'-depsc');
%     print(fullfile(plotdir,[varname '-BodeMag.png']),'-dpng');
%     figure;
%     for k = 1:length(xposPlot)
%         semilogx(freq,unwrap(angle(frTF(indxx(k),:)))*180/pi,'Color',colors(k,:)); 
%         hold on;
%     end
%     for k = 1:length(xposPlot)
%         semilogx(freqSim,unwrap(angle(frSIM(indxx(k),:)))*180/pi,'d','Color',colors(k,:)); 
%     end
%     title(varname);
%     xlabel('Freq');
%     ylabel('Phs');
%     legend(labels1{:},labels2{:},'NumColumns',2,'Location','best', ...
%         'Interpreter','latex');
%     thesisFormat('LineMarkerSize',4);
%     print(fullfile(plotdir,[varname '-BodePhs.eps']),'-depsc');
%     print(fullfile(plotdir,[varname '-BodePhs.png']),'-dpng');
end