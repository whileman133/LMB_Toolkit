% runSensitivityTF_W_kappa.m
%
% Examine the sensitivity of the impedance predicted by the TF model to the
% values W and kappa varied jointly.
%
% -- Changelog --
% 2024.02.03 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath(fullfile('..'));
TB.addpaths;

% Constants.
modelName = 'cellLMO-P2DM';
ff = logspace(-3,5,100);
socPct = 5;
TdegC = 25;
W = 50;
kappap = W/2;  % => Rw(pos) = W/kappap/2 = 1
kappae = W*4;  % => Rw(eff) = W/kappae = 0.25
mult = [1/2; 2];
plotdir = fullfile( ...
    'plots', ...
    sprintf('TFSENS_Warburg_W=%.0f_kappap=%.0fm_%s-%dpct-%ddegC', ...
    W,kappap*1000,modelName,socPct,TdegC));
plotTitle = sprintf( ...
    'sens~W=%.1f~kappae=%.1f~kappap=%.1f', ...
    W,kappae,kappap);
plotLabels = arrayfun( ...
    @(m)sprintf('$%.1f\\bar{W}$, $%.1f\\bar{\\kappa}^\\mathrm{r}$',m,m), ...
    mult,'UniformOutput',false);
plotLabels = [{'Baseline'}; plotLabels(:)];

% Load cell model.
p2dm = loadCellModel(modelName);
wrm = convertCellModel(p2dm,'RLWRM');
lumped = getCellParams(wrm,'TdegC',25);
lumped.const.W = W;
lumped.pos.kappa = kappap;
lumped.eff.kappa = kappae;

% Configure sensitivity study.
sensStudy.defaults = lumped;
sensStudy.joint.multiplier.kappaW.const.W = mult;
sensStudy.joint.multiplier.kappaW.pos.kappa = mult;
sensStudy.joint.multiplier.kappaW.eff.kappa = mult;

% Run sensitivity study.
sensData = fastopt.runSensitivityStudy( ...
    sensStudy,@(params)calcZ(params,ff,socPct,TdegC));


% Plotting ----------------------------------------------------------------

% Make plot directory.
if ~isfolder(plotdir)
    mkdir(plotdir);
end

% Draw plots.
for data = sensData.results
    Z = [data.output.Z];
    Zb = sensData.baseline.Z;
    % subtract out Z(inf) to show how curve shape changes
    Z = Z - Z(end,:);
    Zb = Zb - Zb(end);

    nlabx = 'reZ';
    nlaby = '-ImZ';
    bmlabx = 'freq';
    bmlaby = 'mag';
    bplabx = 'freq';
    bplaby = 'phs';
    ptitle = plotTitle;
    labels = plotLabels;

    figure;
    if strcmp(data.perturbType,'multiplier')
        colororder([0 0 0; winter(size(Z,2))]);
        plot(real(Zb),-imag(Zb),':'); hold on;
        plot(real(Z),-imag(Z));
    else
        colororder(winter(size(Z,2)));
        plot(real(Z),-imag(Z))
    end
    xlabel(nlabx);
    ylabel(nlaby);
    title(ptitle);
    legend(labels,'Location','best','Interpreter','latex');
    setAxesNyquist;
    thesisFormat();
    print(fullfile(plotdir,['Nyq-' data.paramnameEscaped]),'-depsc');
    print(fullfile(plotdir,['Nyq-' data.paramnameEscaped]),'-dpng');

    figure;
    if strcmp(data.perturbType,'multiplier')
        colororder([0 0 0; winter(size(Z,2))]);
        loglog(ff,abs(Zb),':'); hold on;
        loglog(ff,abs(Z));
    else
        colororder(winter(size(Z,2)));
        loglog(ff,abs(Z));
    end
    xlabel(bmlabx);
    ylabel(bmlaby);
    title(ptitle);
    legend(labels,'Location','best','Interpreter','latex');
    xlim([min(ff) max(ff)]);
    thesisFormat;
    print(fullfile(plotdir,['BodeMag-' data.paramnameEscaped]),'-depsc');
    print(fullfile(plotdir,['BodeMag-' data.paramnameEscaped]),'-dpng');

    figure;
    if strcmp(data.perturbType,'multiplier')
        colororder([0 0 0; winter(size(Z,2))]);
        semilogx(ff,unwrap(angle(Zb))*180/pi,':'); hold on;
        semilogx(ff,unwrap(angle(Z))*180/pi);
    else
        colororder(winter(size(Z,2)));
        semilogx(ff,unwrap(angle(Z))*180/pi);
    end
    xlabel(bplabx);
    ylabel(bplaby);
    title(ptitle);
    legend(labels,'Location','best','Interpreter','latex');
    xlim([min(ff) max(ff)]);
    thesisFormat;
    print(fullfile(plotdir,['BodePhs-' data.paramnameEscaped]),'-depsc');
    print(fullfile(plotdir,['BodePhs-' data.paramnameEscaped]),'-dpng');
end

function data = calcZ(params,ff,socPct,TdegC)
    s = 1j*2*pi*ff(:);
    tfdata = tfLMB(s,params, ...
        'socPct',socPct,'TdegC',TdegC,'Calc22',false);
    Z = tfdata.h11.tfZcell();
    Z = Z(:) + params.pkg.R0 + params.pkg.L0*s;
    data.Z = Z;
end