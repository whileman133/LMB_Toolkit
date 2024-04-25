% runProcessGITT.m
%
% Development script for processing GITT data collected by a Gamry
% potentiostat.
%
% -- Changelog --
% 2023.10.10 | Use processGITT.m function | Wes H.
% 2023.09.30 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath('..');
TB.addpaths;

% Constants.
dataFolder = '.';
gittFolder = 'GITT_Cell395530_P25_2m_20m';
dtaName = 'PWRGITT.DTA';
ocpFile = 'OCP_FreshCell_0C01';  % file w/ regressed OCP data
plotDir = 'plots';
outDir = 'processed';
tcPct = 30; % cutoff for finding slope [% relax duration]
gprTheta = linspace(0,1,1000).';
hyp0.sf = 1;
hyp0.ell = 0.05;
hyp0.sn = 0.1;

% Collect data ------------------------------------------------------------
% Load lab OCP.
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle');
ocpData = load(ocpFile);
[~,indTemp] = max(ocpData.study.testTemperatures); % use highest temp.
ocpfit = ocpData.study.tests(indTemp);
ocptest = ocpfit.ocptest;
ocpmodel = MSMR(ocpfit.MSMR);
QtotAh = ocptest.QAh;
warning('on','MATLAB:dispatcher:UnresolvedFunctionHandle');
% !!! Important: We need to compute theta min/max from the voltage limits of 
% the cell, NOT the values in labOCPFit.MSMR, as our lab data may have been
% augmented with additional data collected over wider lithiation range.
ocpTmp = ocpmodel.ocp('voltage',[ocptest.vmin ocptest.vmax]);
theta100 = min(ocpTmp.theta);
theta0 = max(ocpTmp.theta);

% Fetch t, iapp(t), vcell(t).
file = fullfile(dataFolder,gittFolder,dtaName);
gittData = loadGamryDTA(file,'NotesMode','KeyValue');
gittTab = gittData.tables.CURVE;
time = gittTab.Time(:)';
Vcell = gittTab.Vf(:)';
Iapp = -gittTab.Im(:)'; % Gramy uses opposite sign convention
TdegC = mean(gittTab.Temp);
soc0Pct = gittData.notes.soc0Pct;
name = gittData.notes.name;

% Process GITT sequence ---------------------------------------------------
gittData.time = time;
gittData.Iapp = Iapp;
gittData.Vcell = Vcell;
gittData.soc0Pct = soc0Pct;
gittData.TdegC = TdegC;
ocpData.QAh = QtotAh;
ocpData.theta0 = theta0;
ocpData.theta100 = theta100;
ocpData.U0 = ocpmodel.Uj0;
ocpData.omega = ocpmodel.Wj;
ocpData.X = ocpmodel.Xj;
gittOut = processGITT(gittData,ocpData,'tcPct',tcPct);
segments = gittOut.segments;

% Perform GPR -------------------------------------------------------------
gpr = gprBasic( ...
    [segments.theta],log10([segments.DsRel]),gprTheta,hyp0);
est.log10Ds = gpr.optimized.mu;
est.Sigma = gpr.optimized.Sigma;
est.theta = gprTheta;

% Save Results ------------------------------------------------------------
if ~isfolder(outDir)
    mkdir(outDir);
end
save(fullfile(outDir,sprintf('%s.mat',name)),'gittOut','est');

% Do plotting -------------------------------------------------------------
plotFolder = fullfile(plotDir,name);
if ~isfolder(plotFolder)
    mkdir(plotFolder);
end

% GITT sequence.
figure;
plot(time/3600,Vcell);
xlabel('time-hr');
ylabel('vcell');
title('GITT-vcell-v-time');
thesisFormat;
addInset([55 57],[70 3.9]);
addInset([115 120],[10 3.3]);
print(fullfile(plotFolder,'vcell'),'-depsc');
print(fullfile(plotFolder,'vcell'),'-dpng');

% Ds estimates (initial gpr)
% \mathrm{Li}_x\mathrm{Ni}_y\mathrm{Mn}_z\mathrm{Co}_{1-y-z}\mathrm{O}_2
mu = gpr.initial.mu;
s = gpr.initial.Sigma;
figure;
semilogy(gprTheta,10.^mu,'k:'); hold on;
fill([gprTheta;flipud(gprTheta)],10.^[mu+3*sqrt(s);flipud(mu-3*sqrt(s))],...
   'k','EdgeColor','k','FaceAlpha',0.1,'EdgeAlpha',0.3);
semilogy([segments.theta],[segments.DsRel],'ro');
semilogy(gprTheta,10.^mu,'k:');
set(gca,'xdir','reverse');
xlabel('theta');
ylabel('DsRel');
legend('GPR Mean','GPR 3\sigma Bounds','GITT','Location','best');
title('GITT-DsRel-v-theta');
thesisFormat;
print(fullfile(plotFolder,'DsRel-gpr-init'),'-depsc');
print(fullfile(plotFolder,'DsRel-gpr-init'),'-dpng');

% Ds estimates (optimized gpr)
mu = gpr.optimized.mu;
s = gpr.optimized.Sigma;
figure;
semilogy(gprTheta,10.^mu,'k:'); hold on;
fill([gprTheta;flipud(gprTheta)],10.^[mu+3*sqrt(s);flipud(mu-3*sqrt(s))],...
   'k','EdgeColor','k','FaceAlpha',0.1,'EdgeAlpha',0.3);
semilogy([segments.theta],[segments.DsRel],'ro');
semilogy(gprTheta,10.^mu,'k:');
set(gca,'xdir','reverse');
xlabel('theta');
ylabel('DsRel');
legend('GPR Mean','GPR 3\sigma Bounds','GITT','Location','best');
title('GITT-DsRel-v-theta');
thesisFormat;
print(fullfile(plotFolder,'DsRel-gpr-optim'),'-depsc');
print(fullfile(plotFolder,'DsRel-gpr-optim'),'-dpng');

% Raw GITT data.
ind = round(linspace(0,1,21)*(length(segments)-1))+1;
segments = segments(ind);
figure;
l = tiledlayout(7,ceil(length(segments)/7));
l.Title.String = 'GITT';
l.XLabel.String = 'time-hr';
l.YLabel.String = 'vcell';
for k = 1:length(segments)
    seg = segments(k);
    nexttile;
    plot(seg.trelax/3600,seg.Vrelax);
    title(sprintf('%.1f%% SOC',seg.socPct));
end
thesisFormat('FigSizeInches',[7 10]);
print(fullfile(plotFolder,'seg-vcell'),'-depsc');
print(fullfile(plotFolder,'seg-vcell'),'-dpng');
figure;
l = tiledlayout(7,ceil(length(segments)/7));
l.Title.String = 'GITT';
l.XLabel.String = 'sqrt-trelax';
l.YLabel.String = 'vcell';
for k = 1:length(segments)
    seg = segments(k);
    tdiff = sqrt(seg.trelax + seg.tau) - sqrt(seg.trelax);
    nexttile;
    plot(tdiff,seg.Vrelax); hold on;
    l = refline(seg.slope,seg.intercept);
    l.Color = 'k';
    l.LineStyle = ':';
    title(sprintf('%.1f%% SOC',seg.socPct));
end
thesisFormat('FigSizeInches',[7 10]);
print(fullfile(plotFolder,'seg-vcell-sqrt'),'-depsc');
print(fullfile(plotFolder,'seg-vcell-sqrt'),'-dpng');