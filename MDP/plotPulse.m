clear; close all; clc;

% figure;
% load(fullfile('simdata','cellLMO-Lumped-MSMR-modk0-50pct-100mA-qep.mat'));
% idx1 = find(simData.iapp>0,1,'first');
% idx1 = idx1 + 1;
% idx2 = find(simData.time>30,1,'first');
% colors = spring(length(simData.pulse));
% labels1 = arrayfun(@(x)sprintf('%.1fq_{e}^p',x),simData.depletion,'UniformOutput',false);
% for k = 1:length(simData.pulse)
%     plot(simData.time(idx1:idx2), simData.pulse(k).Vcell(idx1:idx2),'-','Color',colors(k,:)); hold on;
% end
% load(fullfile('simdata','cellLMO-Lumped-MSMR-modk0-50pct-100mA-qes.mat'));
% colors = spring(length(simData.pulse));
% labels2 = arrayfun(@(x)sprintf('%.1fq_{e}^s',x),simData.depletion,'UniformOutput',false);
% for k = 1:length(simData.pulse)
%     plot(simData.time(idx1:idx2), simData.pulse(k).Vcell(idx1:idx2),':','Color',colors(k,:)); hold on;
% end
% xlabel('Time [s]');
% ylabel('Cell Voltage [V]');
% title(sprintf('Pulse Response (COMSOL, SOC_0=%d%%, T=%ddegC)',simData.socPct,simData.TdegC));
% legend(labels1{:},labels2{:},'NumColumns',2);
% thesisFormat([0.2 0.1 0.1 0.1]);
% exportgraphics(gca,fullfile('plots','pulse-response-qe.eps'));
% exportgraphics(gca,fullfile('plots','pulse-response-qe.png'));

% figure;
% load(fullfile('simdata','cellLMO-Lumped-MSMR-modk0-50pct-100mA-25degC-psikD.mat'));
% idx1 = find(simData.iapp>0,1,'first');
% idx1 = idx1 + 1;
% idx2 = find(simData.time>30,1,'first');
% colors = spring(length(simData.pulse));
% labels = arrayfun(@(x)sprintf('%.1f(\\psi,\\kappa_D)',x),simData.depletion, ...
%     'UniformOutput',false);
% for k = 1:length(simData.pulse)
%     plot(simData.time(idx1:idx2), simData.pulse(k).Vcell(idx1:idx2),'-', ...
%         'Color',colors(k,:)); hold on;
% end
% load(fullfile('simdata','cellLMO-Lumped-MSMR-modk0-50pct-100mA-25degC-qe-all.mat'));
% colors = spring(length(simData.pulse));
% labels2 = arrayfun(@(x)sprintf('%.1fq_{e}^r',x),simData.depletion,'UniformOutput',false);
% for k = 1:length(simData.pulse)
%     plot(simData.time(idx1:idx2), simData.pulse(k).Vcell(idx1:idx2),'--', ...
%         'Color',colors(k,:)); hold on;
% end
% xlabel('Time [s]');
% ylabel('Cell Voltage [V]');
% title(sprintf('Charge-Nuetral Pulse, Part 1: Discharge (SOC_0=%d%%)', ...
%     simData.socPct));
% legend(labels{:},labels2{:},'NumColumns',2);
% thesisFormat([0.3 0.1 0.3 0.1]);
% exportgraphics(gca,fullfile('plots','pulse-response-psikD.eps'));
% exportgraphics(gca,fullfile('plots','pulse-response-psikD.png'));
% 
% figure;
% load(fullfile('simdata','cellLMO-Lumped-MSMR-modk0-50pct-100mA-25degC-psikD.mat'));
% idx1 = find(simData.iapp<0,1,'first');
% idx1 = idx1 + 1;
% idx2 = find(simData.time>60,1,'first');
% colors = spring(length(simData.pulse));
% labels = arrayfun(@(x)sprintf('%.1f(\\psi,\\kappa_D)',x),simData.depletion, ...
%     'UniformOutput',false);
% for k = 1:length(simData.pulse)
%     plot(simData.time(idx1:idx2), simData.pulse(k).Vcell(idx1:idx2),'-', ...
%         'Color',colors(k,:)); hold on;
% end
% load(fullfile('simdata','cellLMO-Lumped-MSMR-modk0-50pct-100mA-25degC-qe-all.mat'));
% colors = spring(length(simData.pulse));
% labels2 = arrayfun(@(x)sprintf('%.1fq_{e}^r',x),simData.depletion,'UniformOutput',false);
% for k = 1:length(simData.pulse)
%     plot(simData.time(idx1:idx2), simData.pulse(k).Vcell(idx1:idx2),'--', ...
%         'Color',colors(k,:)); hold on;
% end
% xlabel('Time [s]');
% ylabel('Cell Voltage [V]');
% title(sprintf('Charge-Nuetral Pulse, Part 2: Charge (SOC_0=%d%%)', ...
%     simData.socPct));
% legend(labels{:},labels2{:},'NumColumns',2,'Location','best');
% thesisFormat([0.3 0.1 0.3 0.1]);
% exportgraphics(gca,fullfile('plots','pulse-response-psikD2.eps'));
% exportgraphics(gca,fullfile('plots','pulse-response-psikD2.png'));

figure;
load(fullfile('simdata','cellLMO-Lumped-MSMR-modk0-50pct-100mA-25degC-allelec.mat'));
idx1 = find(simData.iapp>0,1,'first');
idx1 = idx1 + 1;
idx2 = find(simData.time>30,1,'first');
colors = spring(length(simData.pulse));
labels = arrayfun(@(x)sprintf('%.1f(\\psi,\\kappa_D,q_{e}^r)',x),simData.depletion, ...
    'UniformOutput',false);
for k = 1:length(simData.pulse)
    plot(simData.time(idx1:idx2), simData.pulse(k).Vcell(idx1:idx2),'-', ...
        'Color',colors(k,:)); hold on;
end
xlabel('Time [s]');
ylabel('Cell Voltage [V]');
title(sprintf('Charge-Nuetral Pulse, Part 1: Discharge (SOC_0=%d%%)', ...
    simData.socPct));
legend(labels{:},'Location','best');
thesisFormat([0.3 0.1 0.3 0.1]);
exportgraphics(gca,fullfile('plots','pulse-response-allelec1.eps'));
exportgraphics(gca,fullfile('plots','pulse-response-allelec1.png'));

figure;
load(fullfile('simdata','cellLMO-Lumped-MSMR-modk0-50pct-100mA-25degC-allelec.mat'));
idx1 = find(simData.iapp<0,1,'first');
idx1 = idx1 + 1;
idx2 = find(simData.time>60,1,'first');
colors = spring(length(simData.pulse));
labels = arrayfun(@(x)sprintf('%.1f(\\psi,\\kappa_D,q_{e}^r)',x),simData.depletion, ...
    'UniformOutput',false);
for k = 1:length(simData.pulse)
    plot(simData.time(idx1:idx2), simData.pulse(k).Vcell(idx1:idx2),'-', ...
        'Color',colors(k,:)); hold on;
end
xlabel('Time [s]');
ylabel('Cell Voltage [V]');
title(sprintf('Charge-Nuetral Pulse, Part 2: Charge (SOC_0=%d%%)', ...
    simData.socPct));
legend(labels{:},'Location','best');
thesisFormat([0.3 0.1 0.3 0.1]);
exportgraphics(gca,fullfile('plots','pulse-response-allelec2.eps'));
exportgraphics(gca,fullfile('plots','pulse-response-allelec2.png'));