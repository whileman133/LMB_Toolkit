% testLayerReductionHypothesis.m
%
% -- Changelog --
% 2023.12.17 | Make code DRY with variable loop
% 2023.04.05 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath(fullfile('..','..'));
TB.addpaths;

% Constants.
cellFile = 'cellLMO-P2DM.xlsx';
plotdir = fullfile('plots','layer-reduction-hypothesis');
freq = logspace(-3.2,5,100);   % Frequency points to evalulate in the spectrum [Hz].
socPct = 5;                    % Cell SOC setpoint [%].
TdegC = 25;                    % Cell temperature [degC].
eEpsDllRel = [0.125 0.25 0.5 1 2];   % Porosity of the dead-lithium layer w/r/t sep porosity.
vars.Thetae = struct('xpos',[0 2 2.5 3]);
%vars.PhieTilde = struct('xpos',[2 2.5 3]); % =0 at x=0, so exclude

% Load cell models.
p2dm = loadCellModel(cellFile);
eEpsSep = getCellParams(p2dm,'sep.eEps');
for k = length(eEpsDllRel):-1:1
    p.dll.eEps = min(1,eEpsSep*eEpsDllRel(k));
    tmp = setCellParam(p2dm,p);
    tmp = convertCellModel(tmp,'WRM');
    modelFull(k) = convertCellModel(tmp,'WRM');
    modelReduced(k) = convertCellModel(tmp,'RLWRM');
end

% Compute frequency response.
for k = length(eEpsDllRel):-1:1
    tfFull(k) = tfLMB(1j*2*pi*freq,modelFull(k), ...
        'socPct',socPct,'TdegC',TdegC);
    tfReduced(k) = tfLMB(1j*2*pi*freq,modelReduced(k), ...
        'socPct',socPct,'TdegC',TdegC);
end

% Genrate directory in which to place plots.
if ~isfolder(plotdir)
    mkdir(plotdir);
end


% Plot linear TFs ---------------------------------------------------------

% Plot linear impedance (Nyquist).
for n = 1:length(tfReduced)
    ZcellFull = tfFull(n).h11.tfZcell();
    ZcellReduced = tfReduced(n).h11.tfZcell();
    magRMSE = rms(abs(ZcellFull)-abs(ZcellReduced));
    phsRMSE = rms(angle(ZcellFull)-angle(ZcellReduced))*180/pi;
    figure;
    plot(real(ZcellFull),-imag(ZcellFull),'m'); hold on;
    plot(real(ZcellReduced),-imag(ZcellReduced),'g:');
    title('Zcell');
    xlabel('Re');
    ylabel('-Im');
    legend('dll + sep','eff only','Location','best');
    annotation('textbox',[0.16, 0.67, 0.175, 0.06], ...
        'String', ...
            ['$\varepsilon_\mathrm{e}^\mathrm{d}' ...
             '/\varepsilon_\mathrm{e}^\mathrm{s}=' ...
             sprintf('%.3f',eEpsDllRel(n)) '$'], ...
        'Interpreter','latex','BackgroundColor','w','FontSize',15);
    annotation('textbox',[0.16, 0.5, 0.175, 0.15], ...
        'String', ...
            [sprintf('RMSE\n') ...
             sprintf('$%.3f\\,\\mathrm{m\\Omega}$ $|Z|$\n',1000*magRMSE) ...
             sprintf('$%.3f^\\circ$ $\\angle Z$',phsRMSE)], ...
        'Interpreter','latex','BackgroundColor','w','FontSize',15);
    setAxesNyquist;
    thesisFormat;
    %ax = addInset([3.451 3.55],[3.25 0.1]);
    %setAxesNyquist('axes',ax,'xdata',[3.451 3.55],'ydata',[0 0.05]);
    print(fullfile(plotdir,sprintf('Zcell-Nyq-eEpsDllRel=%.0f',100*eEpsDllRel(n))),'-depsc');
    print(fullfile(plotdir,sprintf('Zcell-Nyq-eEpsDllRel=%.0f',100*eEpsDllRel(n))),'-dpng');

    % Plot linear TFs for internal variables (Nyqiust).
    varnames = fieldnames(vars);
    for v = 1:length(varnames)
        varname = varnames{v};
        tfname = ['tf' varname];
        vardata = vars.(varname);
        xposPlot = vardata.xpos;
        frFull = tfFull(n).h11.(tfname)(xposPlot);
        frReduced = tfReduced(n).h11.(tfname)(xposPlot);
        labels1 = arrayfun(@(x)sprintf('dll + sep'),xposPlot, ...
            'UniformOutput',false);
        labels2 = arrayfun(@(x)sprintf('eff only $$\\tilde{x}=%.1f$$',x), ...
            xposPlot,'UniformOutput',false);
        colors = cool(length(xposPlot));
    
        figure;
        for k = 1:length(xposPlot)
            plot(real(frFull(:,k)),-imag(frFull(:,k)),'Color',colors(k,:)); 
            hold on;
        end
        for k = 1:length(xposPlot)
            plot(real(frReduced(:,k)),-imag(frReduced(:,k)),':','Color',colors(k,:));
        end
        title(varname);
        xlabel('Re');
        ylabel('-Im');
        legend(labels1{:},labels2{:},'NumColumns',2,'Location','southeast', ...
            'Interpreter','latex');
        annotation('textbox',[0.18, 0.81, 0.175, 0.06], ...
            'String', ...
                ['$\varepsilon_\mathrm{e}^\mathrm{d}' ...
                 '/\varepsilon_\mathrm{e}^\mathrm{s}=' ...
                 sprintf('%.3f',eEpsDllRel(n)) '$'], ...
            'Interpreter','latex','BackgroundColor','w','FontSize',15);
        setAxesNyquist;
        thesisFormat;
        print(fullfile(plotdir,sprintf('%s-Nyq-eEpsDllRel=%.0f',varname,100*eEpsDllRel(n))),'-depsc');
        print(fullfile(plotdir,sprintf('%s-Nyq-eEpsDllRel=%.0f',varname,100*eEpsDllRel(n))),'-dpng');
    end
end