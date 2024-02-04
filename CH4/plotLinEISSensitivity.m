% plotLinEISSensitivity.m

clear; close all; clc;
addpath(fullfile('..'));
TB.addpaths;
modelName = 'cellLMO-P2DM';
p2dm = loadCellModel(modelName);
rlwrm = convertCellModel(p2dm,'RLWRM');

ff = logspace(-3,5,100);
socPct = 5;
TdegC = 25;
lumped = getCellParams(rlwrm,'TdegC',25);
sensStudy.defaults = lumped;
sensStudy.singl.values.pos.alpha = ...
    [0.2 0.8; 0.4 0.8; 0.6 0.8; 0.8 0.8; 0.8 0.6; 0.8 0.4; 0.8 0.2];
sensStudy.singl.multiplier.pos.k0 = [1/5; 1/2; 2; 5];
sensStudy.singl.multiplier.pos.Dsref = [1/5; 1/2; 2; 5];
sensStudy.singl.values.pos.nF = (0.5:0.1:1).';
sensStudy.singl.values.pos.nDL = (0.5:0.1:1).';
sensStudy.singl.multiplier.pos.sigma = [1/100; 1/10; 10; 100];
sensStudy.singl.multiplier.pos.kappa = [1/5; 1/2; 2; 5];
sensStudy.singl.multiplier.sep.kappa = [1/5; 1/2; 2; 5];
sensStudy.singl.multiplier.dll.kappa = [1/5; 1/2; 2; 5];
sensStudy.singl.multiplier.pos.tauW = [0.1; 1/2; 2; 10];
sensStudy.singl.multiplier.sep.tauW = [0.1; 1/2; 2; 10];
sensStudy.singl.multiplier.dll.tauW = [0.1; 1/2; 2; 10];
sensStudy.singl.values.neg.alpha = (0.2:0.2:0.8).';
sensStudy.singl.multiplier.neg.k0 = [1/5; 1/2; 2; 5];
sensStudy.singl.values.neg.nDL = (0.5:0.1:1).';
sensStudy.singl.multiplier.const.psi = [1/5; 1/2; 2; 5];
sensStudy.singl.multiplier.const.W = [1/5; 1/2; 2; 5];
sensStudy.joint.multiplier.psiW.const.psi = [1/5; 1/2; 2; 5];
sensStudy.joint.multiplier.psiW.const.W = [1/5; 1/2; 2; 5];
sensStudy.joint.multiplier.psiWtauW.const.psi = [1/5; 1/2; 2; 5];
sensStudy.joint.multiplier.psiWtauW.const.W = [1/5; 1/2; 2; 5];
sensStudy.joint.multiplier.psiWtauW.pos.tauW = [1/5; 1/2; 2; 5];
sensStudy.joint.multiplier.psiWtauW.sep.tauW = [1/5; 1/2; 2; 5];
sensStudy.joint.multiplier.psiWtauW.dll.tauW = [1/5; 1/2; 2; 5];
sensStudy.joint.multiplier.psitauW.const.psi = [1/5; 1/2; 2; 5];
sensStudy.joint.multiplier.psitauW.pos.tauW = [1/5; 1/2; 2; 5];
sensStudy.joint.multiplier.psitauW.sep.tauW = [1/5; 1/2; 2; 5];
sensStudy.joint.multiplier.psitauW.dll.tauW = [1/5; 1/2; 2; 5];
sensStudy.joint.multiplier.WtauW.const.W = [1/5; 1/2; 2; 5];
sensStudy.joint.multiplier.WtauW.pos.tauW = [1/5; 1/2; 2; 5];
sensStudy.joint.multiplier.WtauW.sep.tauW = [1/5; 1/2; 2; 5];
sensStudy.joint.multiplier.WtauW.dll.tauW = [1/5; 1/2; 2; 5];
sensStudy.joint.multiplier.kappaW.const.W = [1/2; 2; 5; 20];
sensStudy.joint.multiplier.kappaW.pos.kappa = [1/2; 2; 5; 20];
sensStudy.joint.multiplier.kappaW.sep.kappa = [1/2; 2; 5; 20];
sensStudy.joint.multiplier.kappaW.dll.kappa = [1/2; 2; 5; 20];
sensData = fastopt.runSensitivityStudy( ...
    sensStudy,@(params)calcZ(params,ff,socPct,TdegC));

% Make plot directory.
plotdir = fullfile( ...
    'plots', ...
    sprintf('SEN1_%s-%dpct-%ddegC',modelName,socPct,TdegC));
if ~isfolder(plotdir)
    mkdir(plotdir);
end

for data = sensData.results
    Z = [data.output.Z];
    Zb = sensData.baseline.Z;
    if contains(data.basename,{'sigma','kappa'})
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
    if contains(data.basename,{'sigma','kappa'})
        labx = '$(\tilde{Z}_\mathrm{1,1}-\tilde{Z}_\mathrm{1,1}(\infty))''$';
        laby = '$(\tilde{Z}_\mathrm{1,1}-\tilde{Z}_\mathrm{1,1}(\infty))''''$';
    else
        labx = '$\tilde{Z}_\mathrm{1,1}''$';
        laby = '$\tilde{Z}_\mathrm{1,1}''$';
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
    exportgraphics(gcf,fullfile(plotdir,[data.paramnameEscaped '.png']));
    exportgraphics(gcf,fullfile(plotdir,[data.paramnameEscaped '.eps']));
end

function data = calcZ(params,ff,socPct,TdegC)
    tfdata = tfLMB(1j*2*pi*ff,params, ...
        'socPct',socPct,'TdegC',TdegC,'Calc22',false);
    data.Z = tfdata.h11.tfVcell();
end