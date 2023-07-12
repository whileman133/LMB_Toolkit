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
p.const.W = 1;
freq = logspace(-5,5,1000);
s = 1j*2*pi*freq;

tfdata = tfLMB(s,p,'socPct',socPct,'TdegC',TdegC,'Calc22',false);
[~,ZcellData] = tfdata.h11.tfVcell();
Zel = ZcellData.Zel;

L1 = sqrt(p.eff.tauW*s);
Zw1 = +1/p.eff.kappa + 2*(p.const.W/p.eff.kappa)*tanh(L1/2)./L1;
L3 = sqrt(p.pos.tauW*s);
Zw3 = +1/p.pos.kappa/2 + (p.const.W/p.pos.kappa)*tanh(L3/2)./L3;
ZelW = Zw1 + Zw3;

figure;
plot(real(Zel),-imag(Zel),'k'); hold on;
plot(real(ZelW),-imag(ZelW),'r:');
xlabel("Z_{el}' [\Omega]");
ylabel("-Z_{el}'' [\Omega]");
title(sprintf("Impedance contributed by electrolyte (W=%.1f)",p.const.W));
legend('Full-order TF model','FLW approximation');
setAxesNyquist;
thesisFormat;
print(fullfile('plots','APPROX',sprintf('WarburgElectrolyte_W%.0f',p.const.W)),'-depsc');
print(fullfile('plots','APPROX',sprintf('WarburgElectrolyte_W%.0f',p.const.W)),'-dpng');