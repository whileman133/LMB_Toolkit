% plotLinEISSensitivity.m

clear; close all; clc;
addpath(fullfile('..'));
TB.addpaths;
modelName = 'cellLMO-P2DM';
p2dm = loadCellModel(modelName);
wrm = convertCellModel(p2dm,'RLWRM');

ff = logspace(-3,5,100);
socPct = 5;
TdegC = 25;
lumped = getCellParams(wrm,'TdegC',25);
lumped.const.W = 2;
sensStudy.defaults = lumped;
sensStudy.singl.multiplier.pos.kappa = [1/5; 1/2; 2; 5];
sensStudy.singl.multiplier.eff.kappa = [1/5; 1/2; 2; 5];
sensStudy.singl.multiplier.pos.tauW = [0.1; 1/2; 2; 10];
sensStudy.singl.multiplier.eff.tauW = [0.1; 1/2; 2; 10];
sensStudy.singl.multiplier.const.W = [1/5; 1/2; 2; 5];
sensData = fastopt.runSensitivityStudy( ...
    sensStudy,@(params)calcZ(params,ff,socPct,TdegC));

% Make plot directory.
plotdir = fullfile( ...
    'plots', ...
    sprintf('SEN1_EffLayer_%s-%dpct-%ddegC',modelName,socPct,TdegC));
if ~isfolder(plotdir)
    mkdir(plotdir);
end

for data = sensData.results
    Z = [data.output.Z];
    Zb = sensData.baseline.Z;
    subZinf = false; % contains(data.basename,{'sigma','kappa'});
    if subZinf
        % subtract out Z(inf) to show how curve shape changes
        Z = Z - Z(end,:);
        Zb = Zb - Zb(end);
    end

    figure();
    if strcmp(data.perturbType,'multiplier')
        colororder([0 0 0; winter(size(Z,2))]);
        plot(real(Zb),-imag(Zb),':'); hold on;
        plot(real(Z),-imag(Z));
    else
        colororder(winter(size(Z,2)));
        plot(real(Z),-imag(Z))
    end
    if subZinf
        labx = '$(\tilde{Z}_\mathrm{1,1}-\tilde{Z}_\mathrm{1,1}(\infty))''$';
        laby = '$-(\tilde{Z}_\mathrm{1,1}-\tilde{Z}_\mathrm{1,1}(\infty))''''$';
    else
        labx = '$\tilde{Z}_\mathrm{1,1}''$';
        laby = '$-\tilde{Z}_\mathrm{1,1}''$';
    end
    xlabel([labx ...
        ' [$\mathrm{V}\,\mathrm{A}^{-1}$]'],'Interpreter','latex');
    ylabel([laby ...
        ' [$\mathrm{V}\,\mathrm{A}^{-1}$]'],'Interpreter','latex');
    if strcmp(data.analysisType,'joint')
        % Joint paramname too long to include on one line, use abbrev. title!
        title(['$\tilde{v}_\mathrm{cell,1,1}$ ' ...
            'to ' data.paramname],'Interpreter','latex');
    else
        title(['Sensitivity: $\tilde{Z}_\mathrm{1,1}$ ' ...
            'to ' data.paramname ' (Nyquist)'],'Interpreter','latex');
    end
    if strcmp(data.perturbType,'multiplier')
        labels = [{'Baseline'}, data.valuelabels(:)'];
    else
        labels = data.valuelabels; 
    end
    legend(labels,'Location','best','Interpreter','latex');
    setAxesNyquist;
    thesisFormat([0.2 0.1 0.2 0.1]);
    exportgraphics(gcf,fullfile(plotdir,['Nyq-' data.paramnameEscaped '.png']));
    exportgraphics(gcf,fullfile(plotdir,['Nyq-' data.paramnameEscaped '.eps']));

    figure();
    if strcmp(data.perturbType,'multiplier')
        colororder([0 0 0; winter(size(Z,2))]);
        loglog(ff,abs(Zb),':'); hold on;
        loglog(ff,abs(Z));
    else
        colororder(winter(size(Z,2)));
        plot(real(Z),-imag(Z))
    end
    labx = 'Cyclic frequency, $f$ [Hz]';
    if subZinf
        laby = '$|\tilde{Z}_\mathrm{1,1}-\tilde{Z}_\mathrm{1,1}(\infty)|$';
    else
        laby = '$|\tilde{Z}_\mathrm{1,1}|$';
    end
    xlabel(labx,'Interpreter','latex');
    ylabel([laby ...
        ' [$\mathrm{V}\,\mathrm{A}^{-1}$]'],'Interpreter','latex');
    if strcmp(data.analysisType,'joint')
        % Joint paramname too long to include on one line, use abbrev. title!
        title(['$\tilde{v}_\mathrm{cell,1,1}$ ' ...
            'to ' data.paramname],'Interpreter','latex');
    else
        title(['Sensitivity: $\tilde{Z}_\mathrm{1,1}$ ' ...
            'to ' data.paramname ' (Bode Mag.)'],'Interpreter','latex');
    end
    if strcmp(data.perturbType,'multiplier')
        labels = [{'Baseline'}, data.valuelabels(:)'];
    else
        labels = data.valuelabels; 
    end
    legend(labels,'Location','best','Interpreter','latex');
    xlim([min(ff) max(ff)]);
    thesisFormat;
    exportgraphics(gcf,fullfile(plotdir,['BodeMag-' data.paramnameEscaped '.png']));
    exportgraphics(gcf,fullfile(plotdir,['BodeMag-' data.paramnameEscaped '.eps']));

    figure();
    if strcmp(data.perturbType,'multiplier')
        colororder([0 0 0; winter(size(Z,2))]);
        semilogx(ff,unwrap(angle(Zb))*180/pi,':'); hold on;
        semilogx(ff,unwrap(angle(Z))*180/pi);
    else
        colororder(winter(size(Z,2)));
        plot(real(Z),-imag(Z))
    end
    labx = 'Cyclic frequency, $f$ [Hz]';
    if subZinf
        laby = '$\angle(\tilde{Z}_\mathrm{1,1}-\tilde{Z}_\mathrm{1,1}(\infty))$';
    else
        laby = '$\angle\tilde{Z}_\mathrm{1,1}$';
    end
    xlabel(labx,'Interpreter','latex');
    ylabel([laby ...
        ' [deg]'],'Interpreter','latex');
    if strcmp(data.analysisType,'joint')
        % Joint paramname too long to include on one line, use abbrev. title!
        title(['$\tilde{v}_\mathrm{cell,1,1}$ ' ...
            'to ' data.paramname],'Interpreter','latex');
    else
        title(['Sensitivity: $\tilde{Z}_\mathrm{1,1}$ ' ...
            'to ' data.paramname ' (Bode Phs.)'],'Interpreter','latex');
    end
    if strcmp(data.perturbType,'multiplier')
        labels = [{'Baseline'}, data.valuelabels(:)'];
    else
        labels = data.valuelabels; 
    end
    legend(labels,'Location','best','Interpreter','latex');
    xlim([min(ff) max(ff)]);
    thesisFormat;
    exportgraphics(gcf,fullfile(plotdir,['BodePhs-' data.paramnameEscaped '.png']));
    exportgraphics(gcf,fullfile(plotdir,['BodePhs-' data.paramnameEscaped '.eps']));
end

function data = calcZ(params,ff,socPct,TdegC)
    tfdata = tfLMB(1j*2*pi*ff,params, ...
        'socPct',socPct,'TdegC',TdegC,'Calc22',false);
    data.Z = tfdata.h11.tfVcell();
end