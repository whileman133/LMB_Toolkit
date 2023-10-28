% plotDCGain.m

clear; close all; clc;

cellFile = 'cellLMO-P2DM';
socPct = 5;
TdegC = 25;
freq = logspace(-7,6,1000);
s = 1j*2*pi*freq;
pos = linspace(0,3,10);

model = loadCellModel(cellFile);
model = convertCellModel(model,'WRM');
% p.pos.kappa = 1;
% p.sep.kappa = 1;
% p.dll.kappa = 1;
% model = setCellParam(model,p);

tfData = tfLMB(s,model,'socPct',socPct,'TdegC',TdegC);
[Thetae,ThetaeData] = tfData.h11.tfThetae(pos);
[ThetassStar,ThetassStarData] = tfData.h11.tfThetassStar(pos);
[Phie,PhieData] = tfData.h11.tfPhieTilde(pos);
[Ifdl,IfdlData] = tfData.h11.tfIfdl(pos);
[If,IfData] = tfData.h11.tfIf(pos);
[Eta,EtaData] = tfData.h11.tfEta(pos);
[PhiseStar,PhiseStarData] = tfData.h11.tfPhiseStar(pos);
[ZcellStar,ZcellStarData] = tfData.h11.tfZcellStar();

% Plotting ----------------------------------------------------------------
labels = arrayfun(@(x)sprintf('x=%.2f',x),pos,'UniformOutput',false);

figure; colororder(cool(length(pos)));
plot(real(Thetae),-imag(Thetae)); hold on;
plot(ThetaeData.dcGain,0,'k.');
legend(labels{:},'Location','best','NumColumns',2);
xlabel('Re');
ylabel('-Im')
title('$\tilde{\Theta}_\mathrm{e}(\tilde{x},s)/I_\mathrm{app}(s)$','Interpreter','latex');
setAxesNyquist;
thesisFormat;

figure; colororder(cool(length(pos)));
plot(real(Phie),-imag(Phie)); hold on;
plot(PhieData.dcGain,0,'k.');
legend(labels{:},'Location','best','NumColumns',2);
xlabel('Re');
ylabel('-Im')
title('$\tilde{\Phi}_\mathrm{e}(\tilde{x},s)/I_\mathrm{app}(s)$','Interpreter','latex');
setAxesNyquist;
thesisFormat;

figure; colororder(cool(length(pos)));
plot(real(Ifdl),-imag(Ifdl)); hold on;
plot(IfdlData.dcGain,0,'k.');
legend(labels{:},'Location','best','NumColumns',2);
xlabel('Re');
ylabel('-Im')
title('$I_\mathrm{f+dl}(\tilde{x},s)/I_\mathrm{app}(s)$','Interpreter','latex');
setAxesNyquist;
thesisFormat;

figure; colororder(cool(length(pos)));
plot(real(If),-imag(If));  hold on;
plot(IfData.dcGain,0,'k.');
legend(labels{:},'Location','best','NumColumns',2);
xlabel('Re');
ylabel('-Im')
title('$I_\mathrm{f}(\tilde{x},s)/I_\mathrm{app}(s)$','Interpreter','latex');
setAxesNyquist;
thesisFormat;

figure; colororder(cool(length(pos)));
plot(real(Eta),-imag(Eta));  hold on;
plot(EtaData.dcGain,0,'k.');
legend(labels{:},'Location','best','NumColumns',2);
xlabel('Re');
ylabel('-Im')
title('$\eta(\tilde{x},s)/I_\mathrm{app}(s)$','Interpreter','latex');
setAxesNyquist;
thesisFormat;

figure; colororder(cool(length(pos)));
plot(real(PhiseStar),-imag(PhiseStar)); hold on;
plot(PhiseStarData.dcGain,0,'k.');
legend(labels{:},'Location','best','NumColumns',2);
xlabel('Re');
ylabel('-Im')
title('$\tilde{\Phi}_\mathrm{se}^*(\tilde{x},s)/I_\mathrm{app}(s)$','Interpreter','latex');
setAxesNyquist;
thesisFormat;
ax = addInset([-0.24 -0.02],[1 -0.9]);
setAxesNyquist('axes',ax,'xdata',[-0.24 -0.02]);

figure; colororder(cool(length(pos)));
plot(real(ThetassStar),-imag(ThetassStar)); hold on;
plot(ThetassStarData.dcGain,0,'k.');
legend(labels{:},'Location','best','NumColumns',2);
xlabel('Re');
ylabel('-Im')
title('$\Theta_\mathrm{ss}^*(\tilde{x},s)/I_\mathrm{app}(s)$','Interpreter','latex');
setAxesNyquist;
thesisFormat;

figure;
plot(real(ZcellStar),-imag(ZcellStar)); hold on;
plot(ZcellStarData.dcGain,0,'k.');
xlabel('Re');
ylabel('-Im')
title('$Z_\mathrm{cell}^*(s)/I_\mathrm{app}(s)$','Interpreter','latex');
setAxesNyquist;
thesisFormat;