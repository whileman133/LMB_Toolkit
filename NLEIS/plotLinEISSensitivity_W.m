% plotLinEISSensitivity_W.m

clear; close all; clc;
addpath(fullfile('..'));
TB.addpaths;
modelName = 'cellLMO-P2DM';
p2dm = loadCellModel(modelName);
wrm = convertCellModel(p2dm,'WRM');
W = 50;
kappap = W/2;
mult = [1/2; 2;];
plotTitle = sprintf( ...
    'Linear EIS Sensitivity ($\\bar{W}=%.1f$, $\\bar{\\kappa}^\\mathrm{p}=%.1f\\,\\mathrm{\\Omega}^{-1}$)', ...
    W,kappap);
plotLabels = arrayfun( ...
    @(m)sprintf('$%.1f\\bar{W}$, $%.1f\\bar{\\kappa}^\\mathrm{r}$',m,m), ...
    mult,'UniformOutput',false);
plotLabels = [{'Baseline'}; plotLabels(:)];

ff = logspace(-3,5,100);
socPct = 5;
TdegC = 25;
lumped = getCellParams(wrm,'TdegC',25);
lumped.const.W = W;
lumped.pos.kappa = kappap;
sensStudy.defaults = lumped;
sensStudy.joint.multiplier.kappaW.const.W = mult;
sensStudy.joint.multiplier.kappaW.pos.kappa = mult;
sensStudy.joint.multiplier.kappaW.sep.kappa = mult;
sensStudy.joint.multiplier.kappaW.dll.kappa = mult;
sensData = fastopt.runSensitivityStudy( ...
    sensStudy,@(params)calcZ(params,ff,socPct,TdegC));

% Make plot directory.
plotdir = fullfile( ...
    'plots', ...
    sprintf('SEN1_Warburg_W=%.0f_kappap=%.0fm_%s-%dpct-%ddegC', ...
    W,kappap*1000,modelName,socPct,TdegC));
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
        laby = '$-(\tilde{Z}_\mathrm{1,1}-\tilde{Z}_\mathrm{1,1}(\infty))''''$';
    else
        labx = '$\tilde{Z}_\mathrm{1,1}''$';
        laby = '$-\tilde{Z}_\mathrm{1,1}''$';
    end
    xlabel([labx ...
        ' [$\mathrm{V}\,\mathrm{A}^{-1}$]'],'Interpreter','latex');
    ylabel([laby ...
        ' [$\mathrm{V}\,\mathrm{A}^{-1}$]'],'Interpreter','latex');
    title(plotTitle,'Interpreter','latex');
    legend(plotLabels,'Location','best','Interpreter','latex');
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