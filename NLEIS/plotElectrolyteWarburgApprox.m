% plotElectrolyteWarburgApprox.m
%
% Compare FLW approximation for electrolyte impedance to TF model.
%
% -- Changelog --
% 2023.07.11 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath(fullfile('..'));
TB.addpaths;
modelName = 'cellLMO-P2DM';
p2dm = loadCellModel(modelName);
wrm = convertCellModel(p2dm,'RLWRM');

socPct = 60;
TdegC = 25;
p = getCellParams(wrm,'TdegC',25);
p.const.W = 100;
freq = logspace(-6,5,1000);
s = 1j*2*pi*freq;

tfdata = tfLMB(s,p,'socPct',socPct,'TdegC',TdegC,'Calc22',false);
[~,ZcellData] = tfdata.h11.tfVcell();
Zel = ZcellData.Zel;

L1 = sqrt(p.eff.tauW*s);
Zw1 = +1/p.eff.kappa + 2*(p.const.W/p.eff.kappa)*tanh(L1/2)./L1;
L3 = sqrt(p.pos.tauW*s);
Zw3 = +1/p.pos.kappa/2 + (p.const.W/p.pos.kappa)*tanh(L3/2)./L3;
ZelW = Zw1 + Zw3;
f3 = 1/2/pi/p.pos.tauW;
[~,indf3] = min(abs(freq-f3));

figure;
plot(real(Zel),-imag(Zel),'k'); hold on;
plot(real(ZelW),-imag(ZelW),'r:');
xlabel("Z_{el}' [\Omega]");
ylabel("-Z_{el}'' [\Omega]");
title(sprintf("Nyquist: Electrolyte Impedance (W=%.1f)",p.const.W));
legend('Full-order TF model','FLW approximation');
setAxesNyquist;
thesisFormat;
print(fullfile('plots','APPROX',sprintf('WarburgElectrolyte_Nyq_W%.0f',p.const.W)),'-depsc');
print(fullfile('plots','APPROX',sprintf('WarburgElectrolyte_Nyq_W%.0f',p.const.W)),'-dpng');

figure;
loglog(freq,abs(Zel),'k'); hold on;
loglog(freq,abs(ZelW),'r:');
xlabel("Cyclic frequency, f [Hz]");
ylabel("Impedance magnitude, |Z_{el}| [\Omega]");
title(sprintf("Bode Magnitude: Electrolyte Impedance (W=%.1f)",p.const.W));
legend('Full-order TF model','FLW approximation');
thesisFormat;
print(fullfile('plots','APPROX',sprintf('WarburgElectrolyte_BodeMag_W%.0f',p.const.W)),'-depsc');
print(fullfile('plots','APPROX',sprintf('WarburgElectrolyte_BodeMag_W%.0f',p.const.W)),'-dpng');

figure;
semilogx(freq,angle(Zel)*180/pi,'k'); hold on;
semilogx(freq,angle(ZelW)*180/pi,'r:');
xlabel("Cyclic frequency, f [Hz]");
ylabel("Impedance angle, \angleZ_{el} [\circ]");
title(sprintf("Bode Phase: Electrolyte Impedance (W=%.1f)",p.const.W));
legend('Full-order TF model','FLW approximation','Location','best');
thesisFormat;
print(fullfile('plots','APPROX',sprintf('WarburgElectrolyte_BodePhs_W%.0f',p.const.W)),'-depsc');
print(fullfile('plots','APPROX',sprintf('WarburgElectrolyte_BodePhs_W%.0f',p.const.W)),'-dpng');

% figure;
% plot(real(Zw3),-imag(Zw3),'b'); hold on;
% plot(0,0,'k+');
% xlabel("Z_{W}' [\Omega]");
% ylabel("-Z_{W}'' [\Omega]");
% title("Nyquist: Warburg impedance");
% setAxesNyquist;
% thesisFormat;
% print(fullfile('plots','APPROX','FLW-Impedance'),'-depsc');
% print(fullfile('plots','APPROX','FLW-Impedance'),'-dpng');