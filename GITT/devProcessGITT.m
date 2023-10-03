% devProcessGITT.m
%
% Development script for processing GITT data collected by a Gamry
% potentiostat.
%
% -- Changelog --
% 2023.09.30 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath('..');
TB.addpaths;

dataFolder = 'labdata';
gittFolder = 'Sion202309_GITT_Cell395530_P25_30m_180m';
dtaName = 'PWRGITT.DTA';
ocpFile = 'FinalFit-SionFresh_0C01';  % file w/ regressed OCP data
debug = true;
debugOut = fullfile('plots','debug');

% Load lab OCP.
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle');
ocpData = load( ...
    fullfile(TB.const.OCPROOT,'labdata','fitstruct',ocpFile));
[TdegC,indTemp] = max(ocpData.study.testTemperatures); % use highest temp.
ocpfit = ocpData.study.tests(indTemp);
QtotAh = ocpfit.ocptest.QAh;
warning('on','MATLAB:dispatcher:UnresolvedFunctionHandle');

% Fetch t, iapp(t), vcell(t).
file = fullfile(dataFolder,gittFolder,dtaName);
gittData = loadGamryDTA(file,'NotesMode','KeyValue');
gittTab = gittData.tables.CURVE;
time = gittTab.Time(:)';
Vcell = gittTab.Vf(:)';
Iapp = -gittTab.Im(:)'; % Gramy uses opposite sign convention
soc0Pct = gittData.notes.soc0Pct;
name = gittData.notes.name;
debugFolder = fullfile(debugOut,name);
if ~isfolder(debugFolder)
    mkdir(debugFolder);
end

% Segment pulses ----------------------------------------------------------
I = max(abs(Iapp));       % pulse magnitude [A]
relax = abs(Iapp)<I/100;  % logical indicies to relaxation intervals
% For explaination of lines below, let: 
% relax = [ 0  0  0  1  1  1  0  0  0  1  1  1]
% index:    1  2  3  4  5  6  7  8  9 10 11 12
edges = [0 diff(relax)];  % => [ 0  0  0  1  0  0 -1  0  0  1  0  0]
indRelaxStart = find(edges==1);  % => [4 10]
indRelaxEnd = [find(edges==-1)-1 length(edges)];  % => [6 12]
% Create struct array locating the bounds of each segment.
clear segments;
for k = length(indRelaxStart):-1:1
    if k == 1
        indPulseStart = 1; % first pulse starts at time=0
    else
        indPulseStart = indRelaxEnd(k-1)+1;
    end
    indPulseEnd = indRelaxStart(k)-1;
    indPulse = indPulseStart:indPulseEnd;
    indRelax = indRelaxStart(k):indRelaxEnd(k);
    tpulse = time(indPulse);
    tpulse = tpulse - tpulse(1);
    trelax = time(indRelax);
    trelax = trelax - trelax(1);
    tau = tpulse(end) - tpulse(1);
    Iavg = mean(Iapp(indPulse));
    QdisAh = trapz(tpulse,Iapp(indPulse))/3600;
    segments(k).indPulse = indPulse;
    segments(k).indRelax = indRelax;
    segments(k).tau = tau;
    segments(k).trelax = trelax;
    segments(k).Vrelax = Vcell(indRelax);
    segments(k).Iavg = Iavg;
    segments(k).QdisAh = QdisAh;
    segments(k).socPct = [];  % allocate space for later assignment
end
% Determine SOC at start of each relaxation interval.
socPct = soc0Pct;
for k = 1:length(segments)
    socPct = socPct - 100*segments(k).QdisAh/QtotAh;
    segments(k).socPct = socPct;
end
% Debug output.
if debug
    figure;
    plot(time,Vcell); hold on;
    for k = 1:length(segments)
        seg = segments(k);
        xline(time(seg.indPulse(1)),'r');
        xline(time(seg.indPulse(end)),'r--');
        xline(time(seg.indRelax(1)),'k');
        xline(time(seg.indRelax(end)),'k--');
    end
    xlabel('Time [s]');
    ylabel('Cell voltage [V]');
    title(sprintf('%s: Segmentation',name));
    thesisFormat;
    print(fullfile(debugFolder,'vcell'),'-depsc');
    print(fullfile(debugFolder,'vcell'),'-dpng');
end % if debug

% Process segments --------------------------------------------------------
% Debug output.
if debug
    figure;
    l = tiledlayout(7,ceil(length(segments)/7));
    l.Title.String = sprintf('%s',name);
    l.XLabel.Interpreter = 'latex';
    l.XLabel.Interpreter = 'latex';
    l.XLabel.String = 'Time, $t_{relax}$ [s]';
    l.YLabel.Interpreter = 'latex';
    l.YLabel.String = 'Cell voltage, $v_{cell}$ [V]';
    for k = 1:length(segments)
        seg = segments(k);
        nexttile;
        plot(seg.trelax,seg.Vrelax);
        title(sprintf('%.1f%% SOC',seg.socPct));
    end
    thesisFormat('FigSizeInches',[7 10]);
    print(fullfile(debugFolder,'seg-vcell'),'-depsc');
    print(fullfile(debugFolder,'seg-vcell'),'-dpng');
    figure;
    l = tiledlayout(7,ceil(length(segments)/7));
    l.Title.String = sprintf('%s',name);
    l.XLabel.Interpreter = 'latex';
    l.XLabel.String = '$\sqrt{t_{relax}+t_p}-\sqrt{t_{relax}}$ [$\sqrt{s}$]';
    l.YLabel.Interpreter = 'latex';
    l.YLabel.String = 'Cell voltage, $v_{cell}$ [V]';
    for k = 1:length(segments)
        seg = segments(k);
        tdiff = sqrt(seg.trelax + seg.tau) - sqrt(seg.trelax);
        nexttile;
        plot(tdiff,seg.Vrelax);
        title(sprintf('%.1f%% SOC',seg.socPct));
    end
    thesisFormat('FigSizeInches',[7 10]);
    print(fullfile(debugFolder,'seg-vcell-sqrt'),'-depsc');
    print(fullfile(debugFolder,'seg-vcell-sqrt'),'-dpng');
end % if debug