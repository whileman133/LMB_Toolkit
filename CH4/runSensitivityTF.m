% runSensitivityTF.m
%
% Examine the sensitivity of the impedance predicted by the TF model to the
% values of various model parameters.
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
plotdir = fullfile( ...
    'plots', ...
    sprintf('TFSEN_%s-%dpct-%ddegC',modelName,socPct,TdegC));

% Fetch model parameters.
p2dm = loadCellModel(modelName);
rlwrm = convertCellModel(p2dm,'RLWRM');
lumped = getCellParams(rlwrm,'TdegC',25);
lumped.const.W = 1;
lumped.pkg.L0 = 0;

% Configure sensitivity study.
sensStudy.defaults = lumped;
sensStudy.singl.multiplier.const.W = [1/5; 1/2; 2; 5];
sensStudy.singl.multiplier.pos.kappa = [1/5; 1/2; 2; 5];
sensStudy.singl.multiplier.eff.kappa = [1/5; 1/2; 2; 5];
sensStudy.singl.multiplier.pos.tauW = [0.1; 1/2; 2; 10];
sensStudy.singl.multiplier.eff.tauW = [0.1; 1/2; 2; 10];
sensStudy.singl.multiplier.neg.k0 = [1/5; 1/2; 2; 5];
sensStudy.singl.multiplier.pos.k0 = [1/5; 1/2; 2; 5];
sensStudy.singl.values.neg.nDL = (0.5:0.1:1).';
sensStudy.singl.values.pos.nDL = (0.5:0.1:1).';
sensStudy.singl.multiplier.neg.Cdl = [1/5; 1/2; 2; 5];
sensStudy.singl.multiplier.pos.Cdl = [1/5; 1/2; 2; 5];
sensStudy.singl.multiplier.neg.Rdl = [1/100; 1/10; 10; 100];
sensStudy.singl.multiplier.pos.Rdl = [1/5; 1/2; 2; 5];
sensStudy.singl.multiplier.neg.Rf = [1/5; 1/2; 2; 5];
sensStudy.singl.multiplier.pos.Rf = [1/100; 1/10; 10; 100];
sensStudy.singl.multiplier.pos.Dsref = [1/5; 1/2; 2; 5];
sensStudy.singl.values.pos.nF = (0.5:0.1:1).';
sensStudy.singl.multiplier.pos.sigma = [1/100; 1/10; 10; 100];

% Perform sensitivity study.
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
    ptitle = data.paramname;

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
    if strcmp(data.perturbType,'multiplier')
        labels = [{'Baseline'}, data.valuelabels(:)'];
    else
        labels = data.valuelabels; 
    end
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
    if strcmp(data.perturbType,'multiplier')
        labels = [{'Baseline'}, data.valuelabels(:)'];
    else
        labels = data.valuelabels; 
    end
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
    if strcmp(data.perturbType,'multiplier')
        labels = [{'Baseline'}, data.valuelabels(:)'];
    else
        labels = data.valuelabels; 
    end
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