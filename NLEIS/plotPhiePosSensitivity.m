clear; close all; clc;

freq = logspace(-6,3,10000);
params.W = 80;
params.tau = 10;
params.kappa = 10;

sensStudy.defaults = params;
sensStudy.singl.values.W = [1; 5; 10];
sensStudy.singl.values.kappa = [1; 5; 10];
% sensStudy.joint.multiplier.kappaW.W = [1/2; 2; 5; 20];
% sensStudy.joint.multiplier.kappaW.kappa = [1/2; 2; 5; 20];
sensData = fastopt.runSensitivityStudy( ...
    sensStudy,@(params)calcZ(params,freq));


for data = sensData.results
    Z = [data.output.Z];
    Zb = sensData.baseline.Z;
    figure;
    if strcmp(data.perturbType,'multiplier')
        colororder([0 0 0; winter(size(Z,2))]);
        plot(real(Zb),-imag(Zb),':'); hold on;
        plot(real(Z),-imag(Z));
    else
        colororder(winter(size(Z,2)));
        plot(real(Z),-imag(Z))
    end
    labx = '$\tilde{Z}_\mathrm{e,1,1}''$';
    laby = '$-\tilde{Z}_\mathrm{e,1,1}''''$';
    xlabel([labx ...
        ' [$\mathrm{V}\,\mathrm{A}^{-1}$]'],'Interpreter','latex');
    ylabel([laby ...
        ' [$\mathrm{V}\,\mathrm{A}^{-1}$]'],'Interpreter','latex');
    if strcmp(data.analysisType,'joint')
        % Joint paramname too long to include on one line, use abbrev. title!
        title(['$\tilde{v}_\mathrm{e,1,1}$ ' ...
            'to ' data.paramname],'Interpreter','latex');
    else
        title(['Sensitivity: $\tilde{Z}_\mathrm{e,1,1}$ ' ...
            'to ' data.paramname ' (Nyquist)'],'Interpreter','latex');
    end
    if strcmp(data.perturbType,'multiplier')
        labels = [{'Baseline'}, data.valuelabels(:)'];
    else
        labels = data.valuelabels; 
    end
    legend(labels,'Location','best','Interpreter','latex');
    setAxesNyquist;
    thesisFormat;
end


function data = calcZ(p,freq)
    s = 1j*2*pi*freq(:);
    Lambda = sqrt(p.tau*s);
    Phie = -1/2/p.kappa - (p.W/p.kappa)*tanh(Lambda/2)./Lambda;
    data.Z = -Phie;
end