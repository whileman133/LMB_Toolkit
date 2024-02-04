% runSensitivityTF_neg_interface.m
%
% Examine the sensitivity of the impedance predicted by the TF model to the
% values of the negative-electrode surface parameters varied jointly.
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
alpha = [0.2 0.4 0.6 0.8].';
plotdir = fullfile( ...
    'plots', ...
    sprintf('TFSENS__NegInterface_%s-%dpct-%ddegC', ...
    modelName,socPct,TdegC));
plotLabels = arrayfun( ...
    @(m)sprintf('$\\alpha=%.1f$',m), ...
    alpha,'UniformOutput',false);
plotLabels = [{'Baseline'}; plotLabels(:)];

% Load cell model.
p2dm = loadCellModel(modelName);
wrm = convertCellModel(p2dm,'RLWRM');
lumped = getCellParams(wrm,'TdegC',25);
lumped.const.W = 1;

% Configure sensitivity study.
Cdl = lumped.neg.Cdl./alpha;
Rdl = lumped.neg.Rdl.*alpha;
k0 = lumped.neg.k0./alpha;
Rct = TB.const.f(TdegC)./k0;
R0 = (1-alpha).*1./(1./Rct+1./Rdl);
sensStudy.defaults = lumped;
sensStudy.joint.values.negInt.neg.Cdl = Cdl;
sensStudy.joint.values.negInt.neg.Rdl = Rdl;
sensStudy.joint.values.negInt.neg.k0 = k0;
sensStudy.joint.values.negInt.pkg.R0 = R0;

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
    ptitle = data.paramnameEscaped;
    labels = plotLabels;

    figure;
    colororder([0 0 0; winter(size(Z,2))]);
    plot(real(Zb),-imag(Zb),':'); hold on;
    plot(real(Z),-imag(Z));
    xlabel(nlabx);
    ylabel(nlaby);
    title(ptitle);
    legend(labels,'Location','best','Interpreter','latex');
    setAxesNyquist;
    thesisFormat();
    print(fullfile(plotdir,['Nyq-' data.paramnameEscaped]),'-depsc');
    print(fullfile(plotdir,['Nyq-' data.paramnameEscaped]),'-dpng');

    figure;
    colororder([0 0 0; winter(size(Z,2))]);
    loglog(ff,abs(Zb),':'); hold on;
    loglog(ff,abs(Z));
    xlabel(bmlabx);
    ylabel(bmlaby);
    title(ptitle);
    legend(labels,'Location','best','Interpreter','latex');
    xlim([min(ff) max(ff)]);
    thesisFormat;
    print(fullfile(plotdir,['BodeMag-' data.paramnameEscaped]),'-depsc');
    print(fullfile(plotdir,['BodeMag-' data.paramnameEscaped]),'-dpng');

    figure;
    colororder([0 0 0; winter(size(Z,2))]);
    semilogx(ff,unwrap(angle(Zb))*180/pi,':'); hold on;
    semilogx(ff,unwrap(angle(Z))*180/pi);
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