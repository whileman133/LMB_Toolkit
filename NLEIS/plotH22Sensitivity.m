% plotH22Sensitivity.m

clear; close all; clc;
addpath(fullfile("..","TFS"));
addpath(fullfile("..","UTILITY"));
addpath(fullfile("..","XLSX_CELLDEFS"));
addpath(fullfile("..","MAT_CELLDEFS"));
load('cellLMO.mat');

ff = logspace(-3,5,100);
ss = 1j*2*pi*ff;
TdegC = 25;
socPct = 5;

qep = [lumped.pos.qe/100 lumped.pos.qe/10 lumped.pos.qe lumped.pos.qe*10];
Vcell22_qep = zeros(length(ss),length(qep));
for k = 1:length(qep)
    p = lumped;
    p.pos.qe = qep(k);
    data = tfLMB(ss,p, ...
        'socPct',socPct,'TdegC',TdegC,'Calc22',true);
    Vcell22_qep(:,k) = data.h22.tfVcell();
end

beta = [0.3 0.5 0.7];
Vcell22_beta = zeros(length(ss),length(beta));
for k = 1:length(beta)
    p = lumped;
    p.neg.alpha = beta(k);
    data = tfLMB(ss,p, ...
        'socPct',socPct,'TdegC',TdegC,'Calc22',true);
    Vcell22_beta(:,k) = data.h22.tfVcell();
end

betaj = [
    0.2 0.7
    0.2 0.5
    0.5 0.5
    0.7 0.5
    0.7 0.2
];
[nbj,~] = size(betaj);
Vcell22_betaj = zeros(length(ss),nbj);
for k = 1:nbj
    p = lumped;
    p.pos.alpha = betaj(k,:);
    data = tfLMB(ss,p, ...
        'socPct',socPct,'TdegC',TdegC,'Calc22',true);
    Vcell22_betaj(:,k) = data.h22.tfVcell();
end

sigmap = [lumped.pos.sigma/2, lumped.pos.sigma, lumped.pos.sigma*2];
Vcell22_sigmap = zeros(length(ss),length(sigmap));
for k = 1:length(sigmap)
    p = lumped;
    p.pos.sigma = sigmap(k);
    data = tfLMB(ss,p, ...
        'socPct',socPct,'TdegC',TdegC,'Calc22',true);
    Vcell22_sigmap(:,k) = data.h22.tfVcell();
end

kappap = [lumped.pos.kappa/2, lumped.pos.kappa, lumped.pos.kappa*2];
Vcell22_kappap = zeros(length(ss),length(kappap));
for k = 1:length(kappap)
    p = lumped;
    p.pos.kappa = kappap(k);
    data = tfLMB(ss,p, ...
        'socPct',socPct,'TdegC',TdegC,'Calc22',true);
    Vcell22_kappap(:,k) = data.h22.tfVcell();
end

kappas = [lumped.sep.kappa/2, lumped.sep.kappa, lumped.sep.kappa*2];
Vcell22_kappas = zeros(length(ss),length(kappas));
for k = 1:length(kappas)
    p = lumped;
    p.sep.kappa = kappas(k);
    data = tfLMB(ss,p, ...
        'socPct',socPct,'TdegC',TdegC,'Calc22',true);
    Vcell22_kappas(:,k) = data.h22.tfVcell();
end

kappad = [lumped.DL.kappa/2, lumped.DL.kappa, lumped.DL.kappa*2];
Vcell22_kappad = zeros(length(ss),length(kappad));
for k = 1:length(kappad)
    p = lumped;
    p.DL.kappa = kappad(k);
    data = tfLMB(ss,p, ...
        'socPct',socPct,'TdegC',TdegC,'Calc22',true);
    Vcell22_kappad(:,k) = data.h22.tfVcell();
end

% psi = [lumped.const.psi/2, lumped.const.psi, lumped.const.psi*2];
% Vcell22_psi = zeros(length(ss),length(psi));
% for k = 1:length(psi)
%     p = lumped;
%     p.const.psi = psi(k);
%     data = tfLMB(ss,p, ...
%         'socPct',socPct,'TdegC',TdegC,'Calc22',true);
%     Vcell22_psi(:,k) = data.h22.tfVcell();
% end
% 
% kD = [lumped.const.kD/2, lumped.const.kD, lumped.const.kD*2];
% Vcell22_kD = zeros(length(ss),length(kD));
% for k = 1:length(kD)
%     p = lumped;
%     p.const.kD = kD(k);
%     data = tfLMB(ss,p, ...
%         'socPct',socPct,'TdegC',TdegC,'Calc22',true);
%     Vcell22_kD(:,k) = data.h22.tfVcell();
% end
% 
% varp1(1).psi = lumped.const.psi/2;
% varp1(1).kD = lumped.const.kD/2;
% varp1(2).psi = lumped.const.psi;
% varp1(2).kD = lumped.const.kD;
% varp1(3).psi = 2*lumped.const.psi;
% varp1(3).kD = 2*lumped.const.kD;
% Vcell22_psi_kD = zeros(length(ss),length(varp1));
% for k = 1:length(varp1)
%     p = lumped;
%     p.const.psi = varp1(k).psi;
%     p.const.kD = varp1(k).kD;
%     data = tfLMB(ss,p, ...
%         'socPct',socPct,'TdegC',TdegC,'Calc22',true);
%     Vcell22_psi_kD(:,k) = data.h22.tfVcell();
% end
% 
% qep = [lumped.pos.qe/2, lumped.pos.qe, lumped.pos.qe*2];
% Vcell22_qep = zeros(length(ss),length(psi));
% for k = 1:length(psi)
%     p = lumped;
%     p.pos.qe = qep(k);
%     data = tfLMB(ss,p, ...
%         'socPct',socPct,'TdegC',TdegC,'Calc22',true);
%     Vcell22_qep(:,k) = data.h22.tfVcell();
% end
% 
% varp2(1).psi = lumped.const.psi/2;
% varp2(1).qep = lumped.pos.qe/2;
% varp2(2).psi = lumped.const.psi;
% varp2(2).qep = lumped.pos.qe;
% varp2(3).psi = 2*lumped.const.psi;
% varp2(3).qep = 2*lumped.pos.qe;
% Vcell22_psi_qep = zeros(length(ss),length(varp2));
% for k = 1:length(varp2)
%     p = lumped;
%     p.const.psi = varp2(k).psi;
%     p.pos.qe = varp2(k).qep;
%     data = tfLMB(ss,p, ...
%         'socPct',socPct,'TdegC',TdegC,'Calc22',true);
%     Vcell22_psi_qep(:,k) = data.h22.tfVcell();
% end

labels = arrayfun(@(x)sprintf('q_e^p=%.2e',x),qep,'UniformOutput',false);
figure; colororder(cool(length(qep)));
plot(real(Vcell22_qep),-imag(Vcell22_qep));
xlabel('$\tilde{v}_\mathrm{cell,2,2}''$ [$\mathrm{V}\,\mathrm{A}^{-2}$]','Interpreter','latex');
ylabel('-$\tilde{v}_\mathrm{cell,2,2}''''$ [$\mathrm{V}\,\mathrm{A}^{-2}$]','Interpreter','latex');
title('Nyquist: $\tilde{v}_\mathrm{cell,2,2}$','Interpreter','latex');
legend(labels{:},'Location','best');
setAxesNyquist;
thesisFormat([0.1 0.1 0.1 0.1]);
exportgraphics(gcf,'img/sensh22_Vcell_qep.png');
exportgraphics(gcf,'img/sensh22_Vcell_qep.eps');

labels = arrayfun(@(x)sprintf('\\beta^n=%.2f',x),beta,'UniformOutput',false);
figure; colororder(cool(length(beta)));
plot(real(Vcell22_beta),-imag(Vcell22_beta));
xlabel('$\tilde{v}_\mathrm{cell,2,2}''$ [$\mathrm{V}\,\mathrm{A}^{-2}$]','Interpreter','latex');
ylabel('-$\tilde{v}_\mathrm{cell,2,2}''''$ [$\mathrm{V}\,\mathrm{A}^{-2}$]','Interpreter','latex');
title('Nyquist: $\tilde{v}_\mathrm{cell,2,2}$','Interpreter','latex');
legend(labels{:},'Location','best');
setAxesNyquist;
thesisFormat([0.1 0.1 0.1 0.1]);
exportgraphics(gcf,'img/sensh22_Vcell_beta.png');
exportgraphics(gcf,'img/sensh22_Vcell_beta.eps');

labels = cell(1,nbj);
for k = 1:nbj
    labels{k} = sprintf('\\beta^p={%s}',sprintf('%.2f ',betaj(k,:)));
end
figure; colororder(cool(nbj));
plot(real(Vcell22_betaj),-imag(Vcell22_betaj));
xlabel('$\tilde{v}_\mathrm{cell,2,2}''$ [$\mathrm{V}\,\mathrm{A}^{-2}$]','Interpreter','latex');
ylabel('-$\tilde{v}_\mathrm{cell,2,2}''''$ [$\mathrm{V}\,\mathrm{A}^{-2}$]','Interpreter','latex');
title('Nyquist: $\tilde{v}_\mathrm{cell,2,2}$','Interpreter','latex');
legend(labels{:},'Location','best');
setAxesNyquist;
thesisFormat([0.1 0.1 0.1 0.1]);
exportgraphics(gcf,'img/sensh22_Vcell_betaj.png');
exportgraphics(gcf,'img/sensh22_Vcell_betaj.eps');

labels = arrayfun(@(x)sprintf('\\sigma^p=%.2e \\Omega^{-1}',x),sigmap,'UniformOutput',false);
figure; colororder(cool(length(sigmap)));
plot(real(Vcell22_sigmap),-imag(Vcell22_sigmap));
xlabel('$\tilde{v}_\mathrm{cell,2,2}''$ [$\mathrm{V}\,\mathrm{A}^{-2}$]','Interpreter','latex');
ylabel('-$\tilde{v}_\mathrm{cell,2,2}''''$ [$\mathrm{V}\,\mathrm{A}^{-2}$]','Interpreter','latex');
title('Nyquist: $\tilde{v}_\mathrm{cell,2,2}$','Interpreter','latex');
legend(labels{:},'Location','best');
setAxesNyquist;
thesisFormat([0.1 0.1 0.1 0.1]);
exportgraphics(gcf,'img/sensh22_Vcell_sigmap.png');
exportgraphics(gcf,'img/sensh22_Vcell_sigmap.eps');

Z1 = Vcell22_kappap;
labels = arrayfun(@(x)sprintf('\\kappa^p=%.2e \\Omega^{-1}',x),kappap,'UniformOutput',false);
figure; colororder(cool(length(kappap)));
plot(real(Z1),-imag(Z1));
xlabel('$\tilde{v}_\mathrm{cell,2,2}''$ [$\mathrm{V}\,\mathrm{A}^{-2}$]','Interpreter','latex');
ylabel('-$\tilde{v}_\mathrm{cell,2,2}''''$ [$\mathrm{V}\,\mathrm{A}^{-2}$]','Interpreter','latex');
title('Nyquist: $\tilde{v}_\mathrm{cell,2,2}$','Interpreter','latex');
legend(labels{:},'Location','best');
setAxesNyquist;
thesisFormat([0.1 0.1 0.1 0.1]);
exportgraphics(gcf,'img/sensh22_Vcell_kappap.png');
exportgraphics(gcf,'img/sensh22_Vcell_kappap.eps');

Z1 = Vcell22_kappas;
labels = arrayfun(@(x)sprintf('\\kappa^s=%.2e \\Omega^{-1}',x),kappap,'UniformOutput',false);
figure; colororder(cool(length(kappas)));
plot(real(Z1),-imag(Z1));
xlabel('$\tilde{v}_\mathrm{cell,2,2}''$ [$\mathrm{V}\,\mathrm{A}^{-2}$]','Interpreter','latex');
ylabel('-$\tilde{v}_\mathrm{cell,2,2}''''$ [$\mathrm{V}\,\mathrm{A}^{-2}$]','Interpreter','latex');
title('Nyquist: $\tilde{v}_\mathrm{cell,2,2}$','Interpreter','latex');
legend(labels{:},'Location','best');
setAxesNyquist;
thesisFormat([0.1 0.1 0.1 0.1]);
exportgraphics(gcf,'img/sensh22_Vcell_kappas.png');
exportgraphics(gcf,'img/sensh22_Vcell_kappas.eps');

Z1 = Vcell22_kappad;
labels = arrayfun(@(x)sprintf('\\kappa^d=%.2e \\Omega^{-1}',x),kappap,'UniformOutput',false);
figure; colororder(cool(length(kappas)));
plot(real(Z1),-imag(Z1));
xlabel('$\tilde{v}_\mathrm{cell,2,2}''$ [$\mathrm{V}\,\mathrm{A}^{-2}$]','Interpreter','latex');
ylabel('-$\tilde{v}_\mathrm{cell,2,2}''''$ [$\mathrm{V}\,\mathrm{A}^{-2}$]','Interpreter','latex');
title('Nyquist: $\tilde{v}_\mathrm{cell,2,2}$','Interpreter','latex');
legend(labels{:},'Location','best');
setAxesNyquist;
thesisFormat([0.1 0.1 0.1 0.1]);
exportgraphics(gcf,'img/sensh22_Vcell_kappad.png');
exportgraphics(gcf,'img/sensh22_Vcell_kappad.eps');

% labels = arrayfun(@(x)sprintf('\\psi=%.2e V K^{-1}',x),psi,'UniformOutput',false);
% figure; colororder(cool(length(psi)));
% plot(real(Vcell22_psi),-imag(Vcell22_psi));
% xlabel('$\tilde{v}_\mathrm{cell,2,2}''$ [$\mathrm{V}\,\mathrm{A}^{-2}$]','Interpreter','latex');
% ylabel('-$\tilde{v}_\mathrm{cell,2,2}''''$ [$\mathrm{V}\,\mathrm{A}^{-2}$]','Interpreter','latex');
% title('Nyquist: $\tilde{v}_\mathrm{cell,2,2}$','Interpreter','latex');
% legend(labels{:},'Location','best');
% setAxesNyquist;
% thesisFormat([0.1 0.1 0.1 0.1]);
% exportgraphics(gcf,'img/sensh22_Vcell_psi.png');
% exportgraphics(gcf,'img/sensh22_Vcell_psi.eps');
% 
% labels = arrayfun(@(x)sprintf('\\kappa_D=%.2e V K^{-1}',x),kD,'UniformOutput',false);
% figure; colororder(cool(length(kD)));
% plot(real(Vcell22_kD),-imag(Vcell22_kD));
% xlabel('$\tilde{v}_\mathrm{cell,2,2}''$ [$\mathrm{V}\,\mathrm{A}^{-2}$]','Interpreter','latex');
% ylabel('-$\tilde{v}_\mathrm{cell,2,2}''''$ [$\mathrm{V}\,\mathrm{A}^{-2}$]','Interpreter','latex');
% title('Nyquist: $\tilde{v}_\mathrm{cell,2,2}$','Interpreter','latex');
% legend(labels{:},'Location','best');
% setAxesNyquist;
% thesisFormat([0.1 0.1 0.1 0.1]);
% exportgraphics(gcf,'img/sensh22_Vcell_kD.png');
% exportgraphics(gcf,'img/sensh22_Vcell_kD.eps');
% 
% labels = arrayfun(@(x)sprintf('\\psi=%.2e, \\kappa_D=%.2e V K^{-1}',x.psi,x.kD),varp1,'UniformOutput',false);
% figure; colororder(cool(length(varp1)));
% plot(real(Vcell22_psi_kD),-imag(Vcell22_psi_kD));
% xlabel('$\tilde{v}_\mathrm{cell,2,2}''$ [$\mathrm{V}\,\mathrm{A}^{-2}$]','Interpreter','latex');
% ylabel('-$\tilde{v}_\mathrm{cell,2,2}''''$ [$\mathrm{V}\,\mathrm{A}^{-2}$]','Interpreter','latex');
% title('Nyquist: $\tilde{v}_\mathrm{cell,2,2}$','Interpreter','latex');
% legend(labels{:},'Location','best');
% setAxesNyquist;
% thesisFormat([0.1 0.1 0.1 0.1]);
% exportgraphics(gcf,'img/sensh22_Vcell_psi_kD.png');
% exportgraphics(gcf,'img/sensh22_Vcell_psi_kD.eps');
% 
% labels = arrayfun(@(x)sprintf('q_e^p=%.2e Ah',x),qep,'UniformOutput',false);
% figure; colororder(cool(length(qep)));
% plot(real(Vcell22_qep),-imag(Vcell22_qep));
% xlabel('$\tilde{v}_\mathrm{cell,2,2}''$ [$\mathrm{V}\,\mathrm{A}^{-2}$]','Interpreter','latex');
% ylabel('-$\tilde{v}_\mathrm{cell,2,2}''''$ [$\mathrm{V}\,\mathrm{A}^{-2}$]','Interpreter','latex');
% title('Nyquist: $\tilde{v}_\mathrm{cell,2,2}$','Interpreter','latex');
% legend(labels{:},'Location','best');
% setAxesNyquist;
% thesisFormat([0.1 0.1 0.1 0.1]);
% exportgraphics(gcf,'img/sensh22_Vcell_qep.png');
% exportgraphics(gcf,'img/sensh22_Vcell_qep.eps');
% 
% labels = arrayfun(@(x)sprintf('\\psi=%.2e V K^{-1}, q_e^p=%.2e Ah',x.psi,x.qep),varp2,'UniformOutput',false);
% figure; colororder(cool(length(varp2)));
% plot(real(Vcell22_psi_qep),-imag(Vcell22_psi_qep));
% xlabel('$\tilde{v}_\mathrm{cell,2,2}''$ [$\mathrm{V}\,\mathrm{A}^{-2}$]','Interpreter','latex');
% ylabel('-$\tilde{v}_\mathrm{cell,2,2}''''$ [$\mathrm{V}\,\mathrm{A}^{-2}$]','Interpreter','latex');
% title('Nyquist: $\tilde{v}_\mathrm{cell,2,2}$','Interpreter','latex');
% legend(labels{:},'Location','best');
% setAxesNyquist;
% thesisFormat([0.1 0.1 0.1 0.1]);
% exportgraphics(gcf,'img/sensh22_Vcell_psi_qep.png');
% exportgraphics(gcf,'img/sensh22_Vcell_psi_qep.eps');