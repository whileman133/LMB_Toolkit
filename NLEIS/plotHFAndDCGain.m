% plotHFAndDCGain.m
%
% Compare analytic expressions for high-frequency and DC gain to the
% transfer-function model.
%
% -- Changelog --
% 2023.10.31 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;

cellFile = 'cellLMO-P2DM';
socPct = 5;
TdegC = 25;
freq = logspace(-4,6,1000);
s = [0, 1j*2*pi*freq];
pos = linspace(0,3,10);

model = loadCellModel(cellFile);
model = convertCellModel(model,'WRM');

tfData = tfLMB(s,model,'socPct',socPct,'TdegC',TdegC);
[Thetae,ThetaeData] = tfData.h11.tfThetae(pos);
[ThetassStar,ThetassStarData] = tfData.h11.tfThetassStar(pos);
[Phie,PhieData] = tfData.h11.tfPhieTilde(pos);
[Phis,PhisData] = tfData.h11.tfPhisTilde(pos);
[Ifdl,IfdlData] = tfData.h11.tfIfdl(pos);
[If,IfData] = tfData.h11.tfIf(pos);
[Eta,EtaData] = tfData.h11.tfEta(pos);
[PhiseStar,PhiseStarData] = tfData.h11.tfPhiseStar(pos);
[ZcellStar,ZcellStarData] = tfData.h11.tfZcellStar();
Zcell = tfData.h11.tfZcell();

% Plotting ----------------------------------------------------------------
labels = arrayfun(@(x)sprintf('x=%.2f',x),pos,'UniformOutput',false);
markerSize = 3.5;

figure; colororder(cool(length(pos)));
plot(real(Thetae),-imag(Thetae)); hold on;
plot(ThetaeData.dcGain,0,'o');
plot(ThetaeData.hfGain,0,'^');
legend(labels{:},'Location','best','NumColumns',2);
xlabel('Re');
ylabel('-Im')
title('$\tilde{\Theta}_\mathrm{e}(\tilde{x},s)/I_\mathrm{app}(s)$','Interpreter','latex');
setAxesNyquist;
thesisFormat('LineMarkerSize',markerSize);
print(fullfile('plots','dc-hf-gains','Thetae'),'-depsc');
print(fullfile('plots','dc-hf-gains','Thetae'),'-dpng');

figure; colororder(cool(length(pos)));
plot(real(Phie),-imag(Phie)); hold on;
plot(PhieData.dcGain,0,'o');
plot(PhieData.hfGain,0,'^');
legend(labels{:},'Location','best','NumColumns',2);
xlabel('Re');
ylabel('-Im')
title('$\tilde{\Phi}_\mathrm{e}(\tilde{x},s)/I_\mathrm{app}(s)$','Interpreter','latex');
setAxesNyquist;
thesisFormat('LineMarkerSize',markerSize);
print(fullfile('plots','dc-hf-gains','Phie'),'-depsc');
print(fullfile('plots','dc-hf-gains','Phie'),'-dpng');

figure; colororder(cool(length(pos)));
plot(real(Phis),-imag(Phis)); hold on;
plot(PhisData.dcGain,0,'o');
plot(PhisData.hfGain,0,'^');
legend(labels{:},'Location','best','NumColumns',2);
xlabel('Re');
ylabel('-Im')
title('$\tilde{\Phi}_\mathrm{s}(\tilde{x},s)/I_\mathrm{app}(s)$','Interpreter','latex');
%setAxesNyquist;
thesisFormat('LineMarkerSize',markerSize);
print(fullfile('plots','dc-hf-gains','Phis'),'-depsc');
print(fullfile('plots','dc-hf-gains','Phis'),'-dpng');

figure; colororder(cool(length(pos)));
plot(real(Ifdl),-imag(Ifdl)); hold on;
plot(IfdlData.dcGain,0,'o');
plot(IfdlData.hfGain,0,'^');
legend(labels{:},'Location','best','NumColumns',2);
xlabel('Re');
ylabel('-Im')
title('$I_\mathrm{f+dl}(\tilde{x},s)/I_\mathrm{app}(s)$','Interpreter','latex');
setAxesNyquist;
thesisFormat('LineMarkerSize',markerSize);
print(fullfile('plots','dc-hf-gains','Ifdl'),'-depsc');
print(fullfile('plots','dc-hf-gains','Ifdl'),'-dpng');

figure; colororder(cool(length(pos)));
plot(real(If),-imag(If));  hold on;
plot(IfData.dcGain,0,'o');
plot(IfData.hfGain,0,'^');
legend(labels{:},'Location','best','NumColumns',2);
xlabel('Re');
ylabel('-Im')
title('$I_\mathrm{f}(\tilde{x},s)/I_\mathrm{app}(s)$','Interpreter','latex');
setAxesNyquist;
thesisFormat('LineMarkerSize',markerSize);
print(fullfile('plots','dc-hf-gains','If'),'-depsc');
print(fullfile('plots','dc-hf-gains','If'),'-dpng');

figure; colororder(cool(length(pos)));
plot(real(Eta),-imag(Eta));  hold on;
plot(EtaData.dcGain,0,'o');
plot(EtaData.hfGain,0,'^');
legend(labels{:},'Location','best','NumColumns',2);
xlabel('Re');
ylabel('-Im')
title('$\eta(\tilde{x},s)/I_\mathrm{app}(s)$','Interpreter','latex');
setAxesNyquist;
thesisFormat('LineMarkerSize',markerSize);
print(fullfile('plots','dc-hf-gains','Eta'),'-depsc');
print(fullfile('plots','dc-hf-gains','Eta'),'-dpng');

figure; colororder(cool(length(pos)));
plot(real(PhiseStar),-imag(PhiseStar)); hold on;
plot(PhiseStarData.dcGain,0,'o');
plot(PhiseStarData.hfGain,0,'^');
legend(labels{:},'Location','best','NumColumns',2);
xlabel('Re');
ylabel('-Im')
title('$\tilde{\Phi}_\mathrm{se}^*(\tilde{x},s)/I_\mathrm{app}(s)$','Interpreter','latex');
setAxesNyquist;
thesisFormat('LineMarkerSize',markerSize);
ax = addInset([-0.24 -0.02],[1 -0.9]);
setAxesNyquist('axes',ax,'xdata',[-0.24 -0.02]);
print(fullfile('plots','dc-hf-gains','PhiseStar'),'-depsc');
print(fullfile('plots','dc-hf-gains','PhiseStar'),'-dpng');

figure; colororder(cool(length(pos)));
plot(real(ThetassStar),-imag(ThetassStar)); hold on;
plot(ThetassStarData.dcGain,0,'o');
plot(ThetassStarData.hfGain,0,'^');
legend(labels{:},'Location','best','NumColumns',2);
xlabel('Re');
ylabel('-Im')
title('$\Theta_\mathrm{ss}^*(\tilde{x},s)/I_\mathrm{app}(s)$','Interpreter','latex');
setAxesNyquist;
thesisFormat('LineMarkerSize',markerSize);
print(fullfile('plots','dc-hf-gains','ThetassStar'),'-depsc');
print(fullfile('plots','dc-hf-gains','ThetassStar'),'-dpng');

figure;
plot(real(ZcellStar),-imag(ZcellStar)); hold on;
plot(ZcellStarData.dcGain,0,'o');
plot(ZcellStarData.hfGain,0,'^');
xlabel('Re');
ylabel('-Im')
title('$Z_\mathrm{cell}^*(s)/I_\mathrm{app}(s)$','Interpreter','latex');
setAxesNyquist;
thesisFormat('LineMarkerSize',markerSize);
print(fullfile('plots','dc-hf-gains','ZcellStar'),'-depsc');
print(fullfile('plots','dc-hf-gains','ZcellStar'),'-dpng');

figure;
plot(real(Zcell),-imag(Zcell));
xlabel('Re');
ylabel('-Im')
title('$Z_\mathrm{cell}(s)/I_\mathrm{app}(s)$','Interpreter','latex');
setAxesNyquist;
thesisFormat('LineMarkerSize',markerSize);
print(fullfile('plots','dc-hf-gains','Zcell'),'-depsc');
print(fullfile('plots','dc-hf-gains','Zcell'),'-dpng');