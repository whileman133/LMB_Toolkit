% plotNonlinEISSensitivity.m

clear; close all; clc;
addpath(fullfile('..'));
TB.addpaths;
modelName = 'cellLMO-P2DM';
p2dm = loadCellModel(modelName);
wrm = convertCellModel(p2dm,'WRM');

ff = logspace(-3,5,100);
socPct = 5;
TdegC = 25;
lumped = getCellParams(wrm,'TdegC',25);

lumped.pos.k0 = lumped.pos.k0*0.2;
sensStudy.defaults = lumped;
% sensStudy.singl.values.pos.alpha = ...
%     [0.2 0.8; 0.4 0.8; 0.6 0.8; 0.8 0.8; 0.8 0.6; 0.8 0.4; 0.8 0.2];
%sensStudy.singl.values.pos.d2Uocp = [-10; -5; -1; 0; 1; 5; 10];
sensStudy.joint.values.Rct.pos.Rct = 0.0249.*[1/5; 1/2; 1; 2; 5];
sensStudy.joint.values.Rct.pos.Rct2inv = -7./[1/5; 1/2; 1; 2; 5];
% sensStudy.singl.multiplier.pos.k0 = [1/5; 1/2; 2; 5];
% sensStudy.singl.multiplier.pos.Dsref = [1/5; 1/2; 2; 5];
% sensStudy.singl.values.pos.nF = (0.5:0.1:1).';
% sensStudy.singl.values.pos.nDL = (0.5:0.1:1).';
% sensStudy.singl.multiplier.pos.sigma = [1/5; 1/2; 2; 5];
% sensStudy.singl.multiplier.pos.kappa = [1/5; 1/2; 2; 5];
% sensStudy.singl.multiplier.sep.kappa = [1/5; 1/2; 2; 5];
% sensStudy.singl.multiplier.DL.kappa = [1/5; 1/2; 2; 5];
% sensStudy.singl.multiplier.pos.qe = [0.1; 1/2; 2; 10];
% sensStudy.singl.multiplier.sep.qe = [0.1; 1/2; 2; 10];
% sensStudy.singl.multiplier.DL.qe = [0.1; 1/2; 2; 10];
% sensStudy.singl.values.neg.alpha = (0.2:0.2:0.8).';
% sensStudy.singl.multiplier.neg.k0 = [1/5; 1/2; 2; 5];
% sensStudy.singl.values.neg.nDL = (0.5:0.1:1).';
% sensStudy.singl.multiplier.const.psi = [1/5; 1/2; 2; 5];
% sensStudy.singl.multiplier.const.kD = [1/5; 1/2; 2; 5];
% sensStudy.joint.multiplier.psikD.const.psi = [1/5; 1/2; 2; 5];
% sensStudy.joint.multiplier.psikD.const.kD = [1/5; 1/2; 2; 5];
% sensStudy.joint.multiplier.psikDqe.const.psi = [1/5; 1/2; 2; 5];
% sensStudy.joint.multiplier.psikDqe.const.kD = [1/5; 1/2; 2; 5];
% sensStudy.joint.multiplier.psikDqe.pos.qe = [1/5; 1/2; 2; 5];
% sensStudy.joint.multiplier.psikDqe.sep.qe = [1/5; 1/2; 2; 5];
% sensStudy.joint.multiplier.psikDqe.DL.qe = [1/5; 1/2; 2; 5];
sensData = fastopt.runSensitivityStudy( ...
    sensStudy,@(params)calcZ(params,ff,socPct,TdegC));

% Make plot directory.
plotdir = fullfile( ...
    'plots', ...
    sprintf('%s-SENS-H2-%dpct-%ddegC',modelName,socPct,TdegC));
if ~isfolder(plotdir)
    mkdir(plotdir);
end

for data = sensData.results
    Z = [data.output.Z];
    figure;
    if strcmp(data.perturbType,'multiplier')
        colororder([0 0 0; spring(size(Z,2))]);
        plot(real(sensData.baseline.Z),-imag(sensData.baseline.Z),':'); 
        hold on;
        plot(real(Z),-imag(Z));
    else
        colororder(spring(size(Z,2)));
        plot(real(Z),-imag(Z));
    end
    xlabel(['$\tilde{Z}_\mathrm{2,2}''$ ' ...
        '[$\mathrm{V}\,\mathrm{A}^{-2}$]'],'Interpreter','latex');
    ylabel(['-$\tilde{Z}_\mathrm{2,2}''''$ ' ...
        '[$\mathrm{V}\,\mathrm{A}^{-2}$]'],'Interpreter','latex');
    if strcmp(data.analysisType,'joint')
        % Joint paramname too long to include on one line, use abbrev. title!
        title(['$\tilde{Z}_\mathrm{2,2}$ ' ...
            'to ' data.paramname],'Interpreter','latex');
    else
        title(['Sensitivity: $\tilde{Z}_\mathrm{2,2}$ ' ...
            'to ' data.paramname ' (Nyquist)'],'Interpreter','latex');
    end
    if strcmp(data.perturbType,'multiplier')
        labels = [{'Baseline'}, data.valuelabels(:)'];
    else
        labels = data.valuelabels; 
    end
    if strcmp(data.paramnameEscaped,'joint-psikDqe')
        legend(labels,'Location','southwest','Interpreter','latex');
    else
        legend(labels,'Location','best','Interpreter','latex');
    end
    setAxesNyquist;
    thesisFormat([0.2 0.1 0.2 0.1]);
    exportgraphics(gcf,fullfile(plotdir,[data.paramnameEscaped '.png']));
    exportgraphics(gcf,fullfile(plotdir,[data.paramnameEscaped '.eps']));
end

function data = calcZ(params,ff,socPct,TdegC)
    tfdata = tfLMB(1j*2*pi*ff,params, ...
        'socPct',socPct,'TdegC',TdegC,'Calc22',true);
    data.Z = tfdata.h22.tfVcell();
end