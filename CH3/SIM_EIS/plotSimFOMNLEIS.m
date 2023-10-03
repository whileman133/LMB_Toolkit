% plotCOMSOLNonlinEIS.m
%
% Plot results of medium-signal EIS simulation for full-order LMB cell.
%
% -- Changelog --
% 2023.04.05 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath(fullfile('..','..'));
TB.addpaths;

simName = 'cellLMO-P2DM-30mA-socSeries';
seriesData = load(fullfile('simdata',[simName '.mat']));
seriesData = seriesData.simData;

% Constants.
xxThetaePlot = 0:3;
xxThetassPlot = [2 2.5 3];
xxEtaPlot = [2 2.5 3];
fThetaePlot = 0.001;
socPctPlot = 50;

[~,indSOC] = min(abs(seriesData.socPct-socPctPlot));
simData = seriesData.socSeries(indSOC);

% Compute linear and second-harmonic spectra from COMSOL data.
% Also evalulate linear TFs of same variables at same x-locs for 
% comparison to COMSOL simulation.
spectra = processEIS(simData, ...
    'NumHarmonics',2,'EvalLinTF',true,'NumTFFreqPoints',200);
tf = tfLMB(1j*2*pi*spectra.tfFreq,simData.param.cellModel,'Calc22',false, ...
    'TdegC',simData.param.TdegC,'socPct',simData.param.socPct);
xxThetae = spectra.xlocs.Thetae;
xxPhise = spectra.xlocs.Phise;
xxPhie = spectra.xlocs.PhieTilde;
xxThetass = spectra.xlocs.Thetass;
xxEta = spectra.xlocs.Eta;
ZcellSim = [spectra.lin.Zcell];
ZcellTF = [spectra.tf.Zcell];
ThetaeSim = [spectra.lin.Thetae];
ThetaeTF = [spectra.tf.Thetae];
PhiseSim = [spectra.lin.Phise];
PhiseTF = [spectra.tf.Phise];
PhieSim = [spectra.lin.PhieTilde];
PhieTF = [spectra.tf.PhieTilde];
ThetassSim = [spectra.lin.Thetass];
ThetassTF = [spectra.tf.Thetass];
EtaSim = [spectra.lin.Eta];
EtaTF = tf.h11.tfEta(xxEta).';

% Compute second-harmonic spectra.
% Zcell2Sim = [spectra.h2.Zcell];
% Zcell2BVP = tf.h22.tfVcell().';
% Thetae2Sim = [spectra.h2.Thetae];
% Thetae2BVP = tf.h22.tfThetae(xxThetae).';
% Phise2Sim = [spectra.h2.Phise];
% Phise2BVP = tf.h22.tfPhise(xxPhise).';
% Phie2Sim = [spectra.h2.PhieTilde];
% Phie2BVP = tf.h22.tfPhie(xxPhie).';
% Thetass2Sim = [spectra.h2.Thetass];
% Thetass2BVP = tf.h22.tfThetass(xxThetass).';
% Eta2Sim = [spectra.h2.Eta];
% Eta2BVP = tf.h22.tfEta(xxEta).';

% Genrate directory in which to place plots.
plotdir = fullfile('plots',simName);
if ~isfolder(plotdir)
    mkdir(plotdir);
end


% Linear (first harmonic) -------------------------------------------------

% Plot linear impedance (Nyqiust).
figure;
plot(real(ZcellTF),-imag(ZcellTF)); hold on;
plot(real(ZcellSim),-imag(ZcellSim),'d');
title(sprintf('Nyquist: Z_{cell} (%s)',simData.param.cellModel.metadata.cell.name));
xlabel('Z'' [\Omega]');
ylabel('-Z'''' [\Omega]');
legend('TF','COMSOL','Location','best');
setAxesNyquist;
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Zcell-Nyq.eps']));
exportgraphics(gcf,fullfile(plotdir,['Zcell-Nyq.png']));

% Plot linear impedance (Bode).
figure;
loglog(spectra.tfFreq,abs(ZcellTF)); hold on;
loglog(spectra.freq,abs(ZcellSim),'d');
xlim([min(spectra.freq) max(spectra.freq)]);
xlabel('Cyclic Frequency [Hz]');
ylabel('|Z_{cell}| [\Omega]');
title(sprintf('Bode Magnitude: Z_{cell} (%s)',simData.param.cellModel.metadata.cell.name));
legend('TF','COMSOL');
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Zcell-BodeMag.eps']));
exportgraphics(gcf,fullfile(plotdir,['Zcell-BodeMag.png']));
figure;
semilogx(spectra.tfFreq,angle(ZcellTF)*180/pi); hold on;
semilogx(spectra.freq,angle(ZcellSim)*180/pi,'d');
xlim([min(spectra.freq) max(spectra.freq)]);
xlabel('Cyclic Frequency [Hz]');
ylabel('\angleZ_{cell} [deg]');
title(sprintf('Bode Phase: Z_{cell} (%s)',simData.param.cellModel.metadata.cell.name));
legend('TF','COMSOL','Location','best');
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Zcell-BodePhs.eps']));
exportgraphics(gcf,fullfile(plotdir,['Zcell-BodePhs.png']));

% Plot linear Thetae/Iapp TF (Nyqiust and Bode).
[~,indxxThetae] = min(abs(xxThetaePlot(:)'-xxThetae(:)));
labels1 = arrayfun(@(x)sprintf('x=%.1f TF',x),xxThetaePlot, ...
    'UniformOutput',false);
labels2 = arrayfun(@(x)sprintf('x=%.1f FOM',x),xxThetaePlot, ...
    'UniformOutput',false);
colors = cool(length(xxThetaePlot));
figure;
for k = 1:length(xxThetaePlot)
    plot(real(ThetaeTF(indxxThetae(k),:)),-imag(ThetaeTF(indxxThetae(k),:)), ...
        'Color',colors(k,:)); 
    hold on;
end
for k = 1:length(xxThetaePlot)
    plot(real(ThetaeSim(indxxThetae(k),:)),-imag(ThetaeSim(indxxThetae(k),:)),'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
title('Nyquist: \Theta_e(j\omega)/Iapp(j\omega)');
xlabel('Real [A^{-1}]');
ylabel('-Imag [A^{-1}]');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','best');
setAxesNyquist;
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Thetae-Nyq.eps']));
exportgraphics(gcf,fullfile(plotdir,['Thetae-Nyq.png']));
figure;
for k = 1:length(xxThetaePlot)
    loglog(spectra.tfFreq,abs(ThetaeTF(indxxThetae(k),:)), ...
        'Color',colors(k,:));
    hold on;
end
for k = 1:length(xxThetaePlot)
    loglog(spectra.freq,abs(ThetaeSim(indxxThetae(k),:)),'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
xlim([min(spectra.freq) max(spectra.freq)]);
xlabel('Cyclic Frequency, $f$ [Hz]','Interpreter','latex');
ylabel('$|\tilde{\theta}_{e,1,1}|$ [$\mathrm{A}^{-1}$]', ...
    'Interpreter','latex');
title( ...
    sprintf('Bode Magnitude: $\\tilde{\\theta}_{e,1,1}$ (%s)', ...
    simData.param.cellModel.name),'Interpreter','latex');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','best');
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Thetae-BodeMag.eps']));
exportgraphics(gcf,fullfile(plotdir,['Thetae-BodeMag.png']));
figure;
for k = 1:length(xxThetaePlot)
    semilogx(spectra.tfFreq,unwrap(angle(ThetaeTF(indxxThetae(k),:)))*180/pi, ...
        'Color',colors(k,:));
    hold on;
end
for k = 1:length(xxThetaePlot)
    semilogx(spectra.freq,unwrap(angle(ThetaeSim(indxxThetae(k),:)))*180/pi,'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
xlim([min(spectra.freq) max(spectra.freq)]);
xlabel('Cyclic Frequency, $f$ [Hz]','Interpreter','latex');
ylabel('$\angle\tilde{\theta}_{e,1,1}$ [deg]', ...
    'Interpreter','latex');
title( ...
    sprintf('Bode Phase: $\\tilde{\\theta}_{e,1,1}$ (%s)', ...
    simData.param.cellModel.name),'Interpreter','latex');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','best');
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Thetae-BodePhase.eps']));
exportgraphics(gcf,fullfile(plotdir,['Thetae-BodePhase.png']));

return;

% Plot linear Thetae11 phasor vs position.
[~,indfThetae1] = min(abs(fThetaePlot(:)'-spectra.tfFreq(:)));
[~,indfThetae2] = min(abs(fThetaePlot(:)'-spectra.freq(:)));
figure;
plot(xxThetae,real(ThetaeTF(:,indfThetae1)),'r-'); hold on;
plot(xxThetae,imag(ThetaeTF(:,indfThetae1)),'b-');
plot(xxThetae,real(ThetaeSim(:,indfThetae2)),'m--');
plot(xxThetae,imag(ThetaeSim(:,indfThetae2)),'g--');
title(sprintf('$\\tilde{\\theta}_{e,1,1}$ phasor vs. position (f=%.2fmHz)',fThetaePlot*1000), ...
    'Interpreter','latex');
xlabel('Normalized position, $\tilde{x}$ [-]','Interpreter','latex');
ylabel('Phasor value, [$\mathrm{A}^{-1}$]','Interpreter','latex');
legend('$\tilde{\theta}_{e,1,1}''$ TF', ...
    '$\tilde{\theta}_{e,1,1}''''$ TF', ...
    '$\tilde{\theta}_{e,1,1}''$ FOM', ...
    '$\tilde{\theta}_{e,1,1}''''$ FOM', ...
    'NumColumns',2,'Location','best','Interpreter','latex');
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Thetae-Spa.eps']));
exportgraphics(gcf,fullfile(plotdir,['Thetae-Spa.png']));

% Plot linear Phise/Iapp TF (Nyqiust and Bode).
labels1 = arrayfun(@(x)sprintf('x=%.1f TF',x),xxPhise, ...
    'UniformOutput',false);
labels2 = arrayfun(@(x)sprintf('x=%.1f FOM',x),xxPhise, ...
    'UniformOutput',false);
colors = cool(length(xxPhise));
figure;
for k = 1:length(xxPhise)
    plot(real(PhiseTF(k,:)),-imag(PhiseTF(k,:)), ...
        'Color',colors(k,:)); 
    hold on;
end
for k = 1:length(xxPhise)
    plot(real(PhiseSim(k,:)),-imag(PhiseSim(k,:)),'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
title('Nyquist: \Phi_{s,e}(j\omega)/Iapp(j\omega)');
xlabel('Real');
ylabel('-Imag');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','northwest');
setAxesNyquist;
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Phise-Nyq.eps']));
exportgraphics(gcf,fullfile(plotdir,['Phise-Nyq.png']));
figure;
for k = 1:length(xxPhise)
    loglog(spectra.tfFreq,abs(PhiseTF(k,:)), ...
        'Color',colors(k,:));
    hold on;
end
for k = 1:length(xxPhise)
    loglog(spectra.freq,abs(PhiseSim(k,:)),'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
xlim([min(spectra.freq) max(spectra.freq)]);
xlabel('Cyclic Frequency, $f$ [Hz]','Interpreter','latex');
ylabel('$|\tilde{\phi}_{s,e,1,1}|$ [$\mathrm{V}\,\mathrm{A}^{-1}$]', ...
    'Interpreter','latex');
title( ...
    sprintf('Bode Magnitude: $\\tilde{\\phi}_{s,e,1,1}$ (%s)', ...
    simData.param.cellModel.name),'Interpreter','latex');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','best');
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Phise-BodeMag.eps']));
exportgraphics(gcf,fullfile(plotdir,['Phise-BodeMag.png']));
figure;
for k = 1:length(xxPhise)
    semilogx(spectra.tfFreq,unwrap(angle(PhiseTF(k,:)))*180/pi, ...
        'Color',colors(k,:));
    hold on;
end
for k = 1:length(xxPhise)
    semilogx(spectra.freq,unwrap(angle(PhiseSim(k,:)))*180/pi,'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
xlim([min(spectra.freq) max(spectra.freq)]);
xlabel('Cyclic Frequency, $f$ [Hz]','Interpreter','latex');
ylabel('$\angle\tilde{\phi}_{s,e,1,1}$ [deg]', ...
    'Interpreter','latex');
title( ...
    sprintf('Bode Phase: $\\tilde{\\phi}_{s,e,1,1}$ (%s)', ...
    simData.param.cellModel.name),'Interpreter','latex');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','best');
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Phise-BodePhase.eps']));
exportgraphics(gcf,fullfile(plotdir,['Phise-BodePhase.png']));

% Plot linear Phie/Iapp TF (Nyqiust and Bode).
labels1 = arrayfun(@(x)sprintf('x=%.1f TF',x),xxPhie, ...
    'UniformOutput',false);
labels2 = arrayfun(@(x)sprintf('x=%.1f FOM',x),xxPhie, ...
    'UniformOutput',false);
colors = cool(length(xxPhise));
figure;
for k = 1:length(xxPhie)
    plot(real(PhieTF(k,:)),-imag(PhieTF(k,:)), ...
        'Color',colors(k,:)); 
    hold on;
end
for k = 1:length(xxPhie)
    plot(real(PhieSim(k,:)),-imag(PhieSim(k,:)),'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
title('Nyquist: \Phi_{e}(j\omega)/Iapp(j\omega)');
xlabel('Real');
ylabel('-Imag');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','northwest');
setAxesNyquist;
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Phie-Nyq.eps']));
exportgraphics(gcf,fullfile(plotdir,['Phie-Nyq.png']));
figure;
for k = 1:length(xxPhie)
    loglog(spectra.tfFreq,abs(PhieTF(k,:)), ...
        'Color',colors(k,:));
    hold on;
end
for k = 1:length(xxPhie)
    loglog(spectra.freq,abs(PhieSim(k,:)),'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
xlim([min(spectra.freq) max(spectra.freq)]);
xlabel('Cyclic Frequency, $f$ [Hz]','Interpreter','latex');
ylabel('$|\tilde{\phi}_{e,1,1}|$ [$\mathrm{V}\,\mathrm{A}^{-1}$]', ...
    'Interpreter','latex');
title( ...
    sprintf('Bode Magnitude: $\\tilde{\\phi}_{e,1,1}$ (%s)', ...
    simData.param.cellModel.name),'Interpreter','latex');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','best');
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Phie-BodeMag.eps']));
exportgraphics(gcf,fullfile(plotdir,['Phie-BodeMag.png']));
figure;
for k = 1:length(xxPhie)
    semilogx(spectra.tfFreq,unwrap(angle(PhieTF(k,:)))*180/pi, ...
        'Color',colors(k,:));
    hold on;
end
for k = 1:length(xxPhie)
    semilogx(spectra.freq,unwrap(angle(PhieSim(k,:)))*180/pi,'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
xlim([min(spectra.freq) max(spectra.freq)]);
xlabel('Cyclic Frequency, $f$ [Hz]','Interpreter','latex');
ylabel('$\angle\tilde{\phi}_{e,1,1}$ [deg]', ...
    'Interpreter','latex');
title( ...
    sprintf('Bode Phase: $\\tilde{\\phi}_{e,1,1}$ (%s)', ...
    simData.param.cellModel.name),'Interpreter','latex');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','best');
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Phie-BodePhase.eps']));
exportgraphics(gcf,fullfile(plotdir,['Phie-BodePhase.png']));

% Plot linear Thetass/Iapp TF (Nyqiust and Bode).
[~,indxxThetass] = min(abs(xxThetassPlot(:)'-xxThetass(:)));
labels1 = arrayfun(@(x)sprintf('x=%.1f TF',x),xxThetassPlot, ...
    'UniformOutput',false);
labels2 = arrayfun(@(x)sprintf('x=%.1f FOM',x),xxThetassPlot, ...
    'UniformOutput',false);
colors = cool(length(xxThetassPlot));
figure;
for k = 1:length(xxThetassPlot)
    plot(real(ThetassTF(indxxThetass(k),:)),-imag(ThetassTF(indxxThetass(k),:)), ...
        'Color',colors(k,:)); 
    hold on;
end
for k = 1:length(xxThetassPlot)
    plot(real(ThetassSim(indxxThetass(k),:)),-imag(ThetassSim(indxxThetass(k),:)),'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
title('Nyquist: \Theta_{ss}(j\omega)/Iapp(j\omega)');
xlabel('Real [A^{-1}]');
ylabel('-Imag [A^{-1}]');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','northwest');
setAxesNyquist;
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Thetass-Nyq.eps']));
exportgraphics(gcf,fullfile(plotdir,['Thetass-Nyq.png']));
figure;
for k = 1:length(xxThetassPlot)
    loglog(spectra.tfFreq,abs(ThetassTF(indxxThetass(k),:)), ...
        'Color',colors(k,:));
    hold on;
end
for k = 1:length(xxThetassPlot)
    loglog(spectra.freq,abs(ThetassSim(indxxThetass(k),:)),'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
xlim([min(spectra.freq) max(spectra.freq)]);
xlabel('Cyclic Frequency, $f$ [Hz]','Interpreter','latex');
ylabel('$|\tilde{\theta}_{ss,1,1}|$ [$\mathrm{A}^{-1}$]', ...
    'Interpreter','latex');
title( ...
    sprintf('Bode Magnitude: $\\tilde{\\theta}_{ss,1,1}$ (%s)', ...
    simData.param.cellModel.name),'Interpreter','latex');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','best');
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Thetass-BodeMag.eps']));
exportgraphics(gcf,fullfile(plotdir,['Thetass-BodeMag.png']));
figure;
for k = 1:length(xxThetassPlot)
    semilogx(spectra.tfFreq,unwrap(angle(ThetassTF(indxxThetass(k),:)))*180/pi, ...
        'Color',colors(k,:));
    hold on;
end
for k = 1:length(xxThetassPlot)
    semilogx(spectra.freq,unwrap(angle(ThetassSim(indxxThetass(k),:)))*180/pi,'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
xlim([min(spectra.freq) max(spectra.freq)]);
xlabel('Cyclic Frequency, $f$ [Hz]','Interpreter','latex');
ylabel('$\angle\tilde{\theta}_{ss,1,1}$ [deg]', ...
    'Interpreter','latex');
title( ...
    sprintf('Bode Phase: $\\tilde{\\theta}_{ss,1,1}$ (%s)', ...
    simData.param.cellModel.name),'Interpreter','latex');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','best');
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Thetass-BodePhase.eps']));
exportgraphics(gcf,fullfile(plotdir,['Thetass-BodePhase.png']));

% Plot linear Eta/Iapp TF (Nyqiust and Bode).
[~,indxxEta] = min(abs(xxEtaPlot(:)'-xxEta(:)));
labels1 = arrayfun(@(x)sprintf('x=%.1f TF',x),xxEtaPlot, ...
    'UniformOutput',false);
labels2 = arrayfun(@(x)sprintf('x=%.1f FOM',x),xxEtaPlot, ...
    'UniformOutput',false);
colors = cool(length(xxEtaPlot));
figure;
for k = 1:length(xxEtaPlot)
    plot(real(EtaTF(indxxEta(k),:)),-imag(EtaTF(indxxEta(k),:)), ...
        'Color',colors(k,:)); 
    hold on;
end
for k = 1:length(xxEtaPlot)
    plot(real(EtaSim(indxxEta(k),:)),-imag(EtaSim(indxxEta(k),:)),'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
title('Nyquist: \eta(j\omega)/Iapp(j\omega)');
xlabel('Real [V A^{-1}]');
ylabel('-Imag [V A^{-1}]');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','northwest');
setAxesNyquist;
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Eta-Nyq.eps']));
exportgraphics(gcf,fullfile(plotdir,['Eta-Nyq.png']));
figure;
for k = 1:length(xxEtaPlot)
    loglog(spectra.tfFreq,abs(EtaTF(indxxEta(k),:)), ...
        'Color',colors(k,:));
    hold on;
end
for k = 1:length(xxEtaPlot)
    loglog(spectra.freq,abs(EtaSim(indxxEta(k),:)),'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
xlim([min(spectra.freq) max(spectra.freq)]);
xlabel('Cyclic Frequency, $f$ [Hz]','Interpreter','latex');
ylabel('$|\tilde{\eta}_{1,1}|$ [$\mathrm{V}\,\mathrm{A}^{-1}$]', ...
    'Interpreter','latex');
title( ...
    sprintf('Bode Magnitude: $\\tilde{\\eta}_{1,1}$ (%s)', ...
    simData.param.cellModel.name),'Interpreter','latex');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','best');
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Eta-BodeMag.eps']));
exportgraphics(gcf,fullfile(plotdir,['Eta-BodeMag.png']));
figure;
for k = 1:length(xxEtaPlot)
    semilogx(spectra.tfFreq,unwrap(angle(EtaTF(indxxEta(k),:)))*180/pi, ...
        'Color',colors(k,:));
    hold on;
end
for k = 1:length(xxEtaPlot)
    semilogx(spectra.freq,unwrap(angle(EtaSim(indxxEta(k),:)))*180/pi,'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
xlim([min(spectra.freq) max(spectra.freq)]);
xlabel('Cyclic Frequency, $f$ [Hz]','Interpreter','latex');
ylabel('$\angle\tilde{\eta}_{1,1}$ [deg]', ...
    'Interpreter','latex');
title( ...
    sprintf('Bode Phase: $\\tilde{\\eta}_{1,1}$ (%s)', ...
    simData.param.cellModel.name),'Interpreter','latex');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','best');
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Eta-BodePhase.eps']));
exportgraphics(gcf,fullfile(plotdir,['Eta-BodePhase.png']));


% Second harmonic ---------------------------------------------------------

% Plot second-harmonic impedance (Nyqiust).
figure;
plot(real(Zcell2BVP),-imag(Zcell2BVP)); hold on;
plot(real(Zcell2Sim),-imag(Zcell2Sim),'d');
title(sprintf('Nyquist: Z_{cell,2,2} (%s)',simData.param.cellModel.name));
xlabel('Z'' [V A^{-2}]');
ylabel('-Z'''' [V A^{-2}]');
legend('BVP','COMSOL','Location','best');
setAxesNyquist;
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Zcell2-Nyq.eps']));
exportgraphics(gcf,fullfile(plotdir,['Zcell2-Nyq.png']));

% Plot second-harmonic impedance (Bode).
figure;
loglog(spectra.tfFreq,abs(Zcell2BVP)); hold on;
loglog(spectra.freq,abs(Zcell2Sim),'d');
xlim([min(spectra.freq) max(spectra.freq)]);
xlabel('Cyclic Frequency [Hz]');
ylabel('|Z_{cell,2,2}| [V A^{-2}]');
title(sprintf('Bode Magnitude: Z_{cell,2,2} (%s)',simData.param.cellModel.name));
legend('BVP','COMSOL');
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Zcell2-BodeMag.eps']));
exportgraphics(gcf,fullfile(plotdir,['Zcell2-BodeMag.png']));
figure;
semilogx(spectra.tfFreq,unwrap(angle(Zcell2BVP))*180/pi); hold on;
semilogx(spectra.freq,unwrap(angle(Zcell2Sim))*180/pi,'d');
xlim([min(spectra.freq) max(spectra.freq)]);
xlabel('Cyclic Frequency [Hz]');
ylabel('\angleZ_{cell,2,2} [deg]');
title(sprintf('Bode Phase: Z_{cell,2,2} (%s)',simData.param.cellModel.name));
legend('BVP','COMSOL','Location','best');
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Zcell2-BodePhs.eps']));
exportgraphics(gcf,fullfile(plotdir,['Zcell2-BodePhs.png']));

% Plot h2 Thetae (Nyqiust and Bode).
labels1 = arrayfun(@(x)sprintf('x=%.1f BVP',x),xxThetaePlot, ...
    'UniformOutput',false);
labels2 = arrayfun(@(x)sprintf('x=%.1f FOM',x),xxThetaePlot, ...
    'UniformOutput',false);
colors = cool(length(xxThetaePlot));
figure;
for k = 1:length(xxThetaePlot)
    plot(real(Thetae2BVP(indxxThetae(k),:)),-imag(Thetae2BVP(indxxThetae(k),:)), ...
        'Color',colors(k,:)); 
    hold on;
end
for k = 1:length(xxThetaePlot)
    plot(real(Thetae2Sim(indxxThetae(k),:)),-imag(Thetae2Sim(indxxThetae(k),:)),'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
title('Nyquist: $\tilde{\theta}_{\mathrm{e},2,2}$','Interpreter','latex');
xlabel('$\tilde{\theta}_{\mathrm{e},2,2}''$ $[A^{-2}]$','Interpreter','latex');
ylabel('$-\tilde{\theta}_{\mathrm{e},2,2}''''$ $[A^{-2}]$','Interpreter','latex');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','best');
setAxesNyquist;
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Thetae2-Nyq.eps']));
exportgraphics(gcf,fullfile(plotdir,['Thetae2-Nyq.png']));
figure;
for k = 1:length(xxThetaePlot)
    loglog(spectra.tfFreq,abs(Thetae2BVP(indxxThetae(k),:)), ...
        'Color',colors(k,:));
    hold on;
end
for k = 1:length(xxThetaePlot)
    loglog(spectra.freq,abs(Thetae2Sim(indxxThetae(k),:)),'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
ylim([10^(-20) 1]);
xlim([min(spectra.freq) max(spectra.freq)]);
xlabel('Cyclic Frequency, $f$ [Hz]','Interpreter','latex');
ylabel('$|\tilde{\theta}_{e,2,2}|$ [$\mathrm{A}^{-2}$]', ...
    'Interpreter','latex');
title( ...
    sprintf('Bode Magnitude: $\\tilde{\\theta}_{e,2,2}$ (%s)', ...
    simData.param.cellModel.name),'Interpreter','latex');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','best');
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Thetae2-BodeMag.eps']));
exportgraphics(gcf,fullfile(plotdir,['Thetae2-BodeMag.png']));
figure;
for k = 1:length(xxThetaePlot)
    semilogx(spectra.tfFreq,unwrap(angle(Thetae2BVP(indxxThetae(k),:)))*180/pi, ...
        'Color',colors(k,:));
    hold on;
end
for k = 1:length(xxThetaePlot)
    semilogx(spectra.freq,unwrap(angle(Thetae2Sim(indxxThetae(k),:)))*180/pi,'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
ylim([-300 120]);
xlim([min(spectra.freq) max(spectra.freq)]);
xlabel('Cyclic Frequency, $f$ [Hz]','Interpreter','latex');
ylabel('$\angle\tilde{\theta}_{e,2,2}$ [deg]', ...
    'Interpreter','latex');
title( ...
    sprintf('Bode Phase: $\\tilde{\\theta}_{e,2,2}$ (%s)', ...
    simData.param.cellModel.name),'Interpreter','latex');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','best');
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Thetae2-BodePhase.eps']));
exportgraphics(gcf,fullfile(plotdir,['Thetae2-BodePhase.png']));

% Plot second-harmonic Thetae22 phasor vs position.
[~,indfThetae1] = min(abs(fThetaePlot(:)'-spectra.tfFreq(:)));
[~,indfThetae2] = min(abs(fThetaePlot(:)'-spectra.freq(:)));
figure;
plot(xxThetae,real(Thetae2BVP(:,indfThetae1)),'r-'); hold on;
plot(xxThetae,imag(Thetae2BVP(:,indfThetae1)),'b-');
plot(xxThetae,real(Thetae2Sim(:,indfThetae2)),'m--');
plot(xxThetae,imag(Thetae2Sim(:,indfThetae2)),'g--');
title(sprintf('$\\tilde{\\theta}_{e,2,2}$ phasor vs. position (f=%.2fmHz)',fThetaePlot*1000), ...
    'Interpreter','latex');
xlabel('Normalized position, $\tilde{x}$ [-]','Interpreter','latex');
ylabel('Phasor value, [$\mathrm{A}^{-2}$]','Interpreter','latex');
legend('$\tilde{\theta}_{e,2,2}''$ BVP', ...
    '$\tilde{\theta}_{e,2,2}''''$ BVP', ...
    '$\tilde{\theta}_{e,2,2}''$ FOM', ...
    '$\tilde{\theta}_{e,2,2}''''$ FOM', ...
    'NumColumns',2,'Location','best','Interpreter','latex');
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Thetae2-Spa.eps']));
exportgraphics(gcf,fullfile(plotdir,['Thetae2-Spa.png']));

% Plot second-harmonic Phise22 phasor (Nyqiust and Bode).
labels1 = arrayfun(@(x)sprintf('x=%.1f BVP',x),xxPhise, ...
    'UniformOutput',false);
labels2 = arrayfun(@(x)sprintf('x=%.1f FOM',x),xxPhise, ...
    'UniformOutput',false);
colors = cool(length(xxPhise));
figure;
for k = 1:length(xxPhise)
    plot(real(Phise2BVP(k,:)),-imag(Phise2BVP(k,:)), ...
        'Color',colors(k,:)); 
    hold on;
end
for k = 1:length(xxPhise)
    plot(real(Phise2Sim(k,:)),-imag(Phise2Sim(k,:)),'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
title('Nyquist: $\tilde{\phi}_{\mathrm{s,e},2,2}$ Phasor','Interpreter','latex');
xlabel('$\tilde{\phi}_{\mathrm{s,e},2,2}''$ $[\mathrm{V}\,\mathrm{A}^{-2}]$','Interpreter','latex');
ylabel('$-\tilde{\phi}_{\mathrm{s,e},2,2}''''$ $[\mathrm{V}\,\mathrm{A}^{-2}]$','Interpreter','latex');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','best');
setAxesNyquist;
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Phise2-Nyq.eps']));
exportgraphics(gcf,fullfile(plotdir,['Phise2-Nyq.png']));
figure;
for k = 1:length(xxPhise)
    loglog(spectra.tfFreq,abs(Phise2BVP(k,:)), ...
        'Color',colors(k,:));
    hold on;
end
for k = 1:length(xxPhise)
    loglog(spectra.freq,abs(Phise2Sim(k,:)),'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
xlim([min(spectra.freq) max(spectra.freq)]);
xlabel('Cyclic Frequency, $f$ [Hz]','Interpreter','latex');
ylabel('$|\tilde{\phi}_{s,e,2,2}|$ [$\mathrm{V}\,\mathrm{A}^{-2}$]', ...
    'Interpreter','latex');
title( ...
    sprintf('Bode Magnitude: $\\tilde{\\phi}_{s,e,2,2}$ (%s)', ...
    simData.param.cellModel.name),'Interpreter','latex');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','best');
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Phise2-BodeMag.eps']));
exportgraphics(gcf,fullfile(plotdir,['Phise2-BodeMag.png']));
figure;
for k = 1:length(xxPhise)
    semilogx(spectra.tfFreq,unwrap(angle(Phise2BVP(k,:)))*180/pi, ...
        'Color',colors(k,:));
    hold on;
end
for k = 1:length(xxPhise)
    semilogx(spectra.freq,unwrap(angle(Phise2Sim(k,:)))*180/pi,'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
xlim([min(spectra.freq) max(spectra.freq)]);
xlabel('Cyclic Frequency, $f$ [Hz]','Interpreter','latex');
ylabel('$\angle\tilde{\phi}_{s,e,2,2}$ [deg]', ...
    'Interpreter','latex');
title( ...
    sprintf('Bode Phase: $\\tilde{\\phi}_{s,e,2,2}$ (%s)', ...
    simData.param.cellModel.name),'Interpreter','latex');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','best');
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Phise2-BodePhase.eps']));
exportgraphics(gcf,fullfile(plotdir,['Phise2-BodePhase.png']));

% Plot second-harmonic Phie22 phasor (Nyqiust and Bode).
labels1 = arrayfun(@(x)sprintf('x=%.1f BVP',x),xxPhie, ...
    'UniformOutput',false);
labels2 = arrayfun(@(x)sprintf('x=%.1f FOM',x),xxPhie, ...
    'UniformOutput',false);
colors = cool(length(xxPhise));
figure;
for k = 1:length(xxPhie)
    plot(real(Phie2BVP(k,:)),-imag(Phie2BVP(k,:)), ...
        'Color',colors(k,:)); 
    hold on;
end
for k = 1:length(xxPhie)
    plot(real(Phie2Sim(k,:)),-imag(Phie2Sim(k,:)),'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
title('Nyquist: $\tilde{\phi}_{\mathrm{e},2,2}$ Phasor','Interpreter','latex');
xlabel('$\tilde{\phi}_{\mathrm{e},2,2}''$ $[\mathrm{V}\,\mathrm{A}^{-2}]$','Interpreter','latex');
ylabel('$-\tilde{\phi}_{\mathrm{e},2,2}''''$ $[\mathrm{V}\,\mathrm{A}^{-2}]$','Interpreter','latex');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','best');
setAxesNyquist;
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Phie2-Nyq.eps']));
exportgraphics(gcf,fullfile(plotdir,['Phie2-Nyq.png']));
figure;
for k = 1:length(xxPhie)
    loglog(spectra.tfFreq,abs(Phie2BVP(k,:)), ...
        'Color',colors(k,:));
    hold on;
end
for k = 1:length(xxPhie)
    loglog(spectra.freq,abs(Phie2Sim(k,:)),'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
xlim([min(spectra.freq) max(spectra.freq)]);
xlabel('Cyclic Frequency, $f$ [Hz]','Interpreter','latex');
ylabel('$|\tilde{\phi}_{e,2,2}|$ [$\mathrm{V}\,\mathrm{A}^{-2}$]', ...
    'Interpreter','latex');
title( ...
    sprintf('Bode Magnitude: $\\tilde{\\phi}_{e,2,2}$ (%s)', ...
    simData.param.cellModel.name),'Interpreter','latex');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','best');
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Phie2-BodeMag.eps']));
exportgraphics(gcf,fullfile(plotdir,['Phie2-BodeMag.png']));
figure;
for k = 1:length(xxPhie)
    semilogx(spectra.tfFreq,unwrap(angle(Phie2BVP(k,:)))*180/pi, ...
        'Color',colors(k,:));
    hold on;
end
for k = 1:length(xxPhie)
    semilogx(spectra.freq,unwrap(angle(Phie2Sim(k,:)))*180/pi,'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
xlim([min(spectra.freq) max(spectra.freq)]);
xlabel('Cyclic Frequency, $f$ [Hz]','Interpreter','latex');
ylabel('$\angle\tilde{\phi}_{e,2,2}$ [deg]', ...
    'Interpreter','latex');
title( ...
    sprintf('Bode Phase: $\\tilde{\\phi}_{e,2,2}$ (%s)', ...
    simData.param.cellModel.name),'Interpreter','latex');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','best');
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Phie2-BodePhase.eps']));
exportgraphics(gcf,fullfile(plotdir,['Phie2-BodePhase.png']));

% Plot h2 Thetass (Nyqiust and Bode).
[~,indxxThetass] = min(abs(xxThetassPlot(:)'-xxThetass(:)));
labels1 = arrayfun(@(x)sprintf('x=%.1f TF',x),xxThetassPlot, ...
    'UniformOutput',false);
labels2 = arrayfun(@(x)sprintf('x=%.1f FOM',x),xxThetassPlot, ...
    'UniformOutput',false);
colors = cool(length(xxThetassPlot));
figure;
for k = 1:length(xxThetassPlot)
    plot(real(Thetass2BVP(indxxThetass(k),:)),-imag(Thetass2BVP(indxxThetass(k),:)), ...
        'Color',colors(k,:)); 
    hold on;
end
for k = 1:length(xxThetassPlot)
    plot(real(Thetass2Sim(indxxThetass(k),:)),-imag(Thetass2Sim(indxxThetass(k),:)),'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
title( ...
    sprintf('Nyquist: $\\tilde{\\theta}_{ss,2,2}$ (%s)', ...
    simData.param.cellModel.name),'Interpreter','latex');
xlabel('$\tilde{\theta}_{\mathrm{ss},2,2}''$ $[\mathrm{A}^{-2}]$','Interpreter','latex');
ylabel('$-\tilde{\theta}_{\mathrm{ss},2,2}''''$ $[\mathrm{A}^{-2}]$','Interpreter','latex');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','best');
setAxesNyquist;
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Thetass2-Nyq.eps']));
exportgraphics(gcf,fullfile(plotdir,['Thetass2-Nyq.png']));
figure;
for k = 1:length(xxThetassPlot)
    loglog(spectra.tfFreq,abs(Thetass2BVP(indxxThetass(k),:)), ...
        'Color',colors(k,:));
    hold on;
end
for k = 1:length(xxThetassPlot)
    loglog(spectra.freq,abs(Thetass2Sim(indxxThetass(k),:)),'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
xlim([min(spectra.freq) max(spectra.freq)]);
xlabel('Cyclic Frequency, $f$ [Hz]','Interpreter','latex');
ylabel('$|\tilde{\theta}_{ss,2,2}|$ [$\mathrm{A}^{-2}$]', ...
    'Interpreter','latex');
title( ...
    sprintf('Bode Magnitude: $\\tilde{\\theta}_{ss,2,2}$ (%s)', ...
    simData.param.cellModel.name),'Interpreter','latex');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','best');
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Thetass2-BodeMag.eps']));
exportgraphics(gcf,fullfile(plotdir,['Thetass2-BodeMag.png']));
figure;
for k = 1:length(xxThetassPlot)
    semilogx(spectra.tfFreq,unwrap(angle(Thetass2BVP(indxxThetass(k),:)))*180/pi, ...
        'Color',colors(k,:));
    hold on;
end
for k = 1:length(xxThetassPlot)
    semilogx(spectra.freq,unwrap(angle(Thetass2Sim(indxxThetass(k),:)))*180/pi,'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
xlim([min(spectra.freq) max(spectra.freq)]);
xlabel('Cyclic Frequency, $f$ [Hz]','Interpreter','latex');
ylabel('$\angle\tilde{\theta}_{ss,2,2}$ [deg]', ...
    'Interpreter','latex');
title( ...
    sprintf('Bode Phase: $\\tilde{\\theta}_{ss,2,2}$ (%s)', ...
    simData.param.cellModel.name),'Interpreter','latex');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','best');
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Thetass2-BodePhase.eps']));
exportgraphics(gcf,fullfile(plotdir,['Thetass2-BodePhase.png']));

% Plot h2 Eta (Nyqiust and Bode).
[~,indxxEta] = min(abs(xxEtaPlot(:)'-xxEta(:)));
labels1 = arrayfun(@(x)sprintf('x=%.1f TF',x),xxEtaPlot, ...
    'UniformOutput',false);
labels2 = arrayfun(@(x)sprintf('x=%.1f FOM',x),xxEtaPlot, ...
    'UniformOutput',false);
colors = cool(length(xxEtaPlot));
figure;
for k = 1:length(xxEtaPlot)
    plot(real(Eta2BVP(indxxEta(k),:)),-imag(Eta2BVP(indxxEta(k),:)), ...
        'Color',colors(k,:)); 
    hold on;
end
for k = 1:length(xxEtaPlot)
    plot(real(Eta2Sim(indxxEta(k),:)),-imag(Eta2Sim(indxxEta(k),:)),'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
title( ...
    sprintf('Nyquist: $\\tilde{\\eta}_{2,2}$ (%s)', ...
    simData.param.cellModel.name),'Interpreter','latex');
xlabel('$\tilde{\eta}_{2,2}''$ $[\mathrm{V}\,\mathrm{A}^{-2}]$','Interpreter','latex');
ylabel('$-\tilde{\eta}_{2,2}''''$ $[\mathrm{V}\,\mathrm{A}^{-2}]$','Interpreter','latex');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','best');
setAxesNyquist;
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Eta2-Nyq.eps']));
exportgraphics(gcf,fullfile(plotdir,['Eta2-Nyq.png']));
figure;
for k = 1:length(xxEtaPlot)
    loglog(spectra.tfFreq,abs(Eta2BVP(indxxEta(k),:)), ...
        'Color',colors(k,:));
    hold on;
end
for k = 1:length(xxEtaPlot)
    loglog(spectra.freq,abs(Eta2Sim(indxxEta(k),:)),'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
xlim([min(spectra.freq) max(spectra.freq)]);
xlabel('Cyclic Frequency, $f$ [Hz]','Interpreter','latex');
ylabel('$|\tilde{\eta}_{2,2}|$ [$\mathrm{V}\,\mathrm{A}^{-1}$]', ...
    'Interpreter','latex');
title( ...
    sprintf('Bode Magnitude: $\\tilde{\\eta}_{2,2}$ (%s)', ...
    simData.param.cellModel.name),'Interpreter','latex');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','best');
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Eta2-BodeMag.eps']));
exportgraphics(gcf,fullfile(plotdir,['Eta2-BodeMag.png']));
figure;
for k = 1:length(xxEtaPlot)
    semilogx(spectra.tfFreq,unwrap(angle(Eta2BVP(indxxEta(k),:)))*180/pi, ...
        'Color',colors(k,:));
    hold on;
end
for k = 1:length(xxEtaPlot)
    semilogx(spectra.freq,unwrap(angle(Eta2Sim(indxxEta(k),:)))*180/pi,'d', ...
        'Color',colors(k,:),'MarkerFaceColor','auto');
end
xlim([min(spectra.freq) max(spectra.freq)]);
xlabel('Cyclic Frequency, $f$ [Hz]','Interpreter','latex');
ylabel('$\angle\tilde{\eta}_{2,2}$ [deg]', ...
    'Interpreter','latex');
title( ...
    sprintf('Bode Phase: $\\tilde{\\eta}_{2,2}$ (%s)', ...
    simData.param.cellModel.name),'Interpreter','latex');
legend(labels1{:},labels2{:},'NumColumns',2,'Location','best');
thesisFormat([0.2 0.1 0.1 0.1]);
exportgraphics(gcf,fullfile(plotdir,['Eta2-BodePhase.eps']));
exportgraphics(gcf,fullfile(plotdir,['Eta2-BodePhase.png']));