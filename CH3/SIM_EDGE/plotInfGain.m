% plotInfGain.m
%
% Compare analytic expressions for infinute-frequency gain to the
% transfer-function model.
%
% -- Changelog --
% 2024.01.06 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;

% Constants.
cellFile = 'cellLMO-P2DM';
socPct = 10;
TdegC = 25;
freq = logspace(-4,6,1000);
s = [0, 1j*2*pi*freq];
xpos = linspace(0,3,10);
vars.Thetae = struct('xpos',[0 1 2 2.5 3]);
vars.ThetassStar = struct('xpos',[2 2.5 3]);
vars.PhieTilde = struct('xpos',[1 2 2.5 3]); % =0 at x=0, so exlcude
vars.PhisTilde = struct('xpos',[2 2.5]); % =0 at x=3, so exclude
vars.PhiseStar = struct('xpos',[2 2.5 3]);
vars.Eta = struct('xpos',[2 2.5 3]);
vars.If = struct('xpos',[2 2.5 3]);
vars.Ifdl = struct('xpos',[2 2.5 3]);
plotdir = fullfile('plots','inf-gain',cellFile);

% Fetch cell model.
model = loadCellModel(cellFile);
model = convertCellModel(model,'WRM');

tfData = tfLMB(s,model,'socPct',socPct,'TdegC',TdegC);
[ZcellStar,ZcellStarData] = tfData.h11.tfZcellStar();

% Plotting ----------------------------------------------------------------
markerSize = 3.5;
if ~isfolder(plotdir)
    mkdir(plotdir);
end

figure;
plot(real(ZcellStar),-imag(ZcellStar)); hold on;
plot(ZcellStarData.dcGain,0,'o');
plot(ZcellStarData.hfGain,0,'^');
xlabel('Re');
ylabel('-Im')
title('ZcellStar');
setAxesNyquist;
thesisFormat('LineMarkerSize',markerSize);
print(fullfile(plotdir,'ZcellStar'),'-depsc');
print(fullfile(plotdir,'ZcellStar'),'-dpng');

varnames = fieldnames(vars);
for k = 1:length(varnames)
    varname = varnames{k};
    vardata = vars.(varname);
    xpos = vardata.xpos;
    [fr, data] = tfData.h11.(['tf' varname])(xpos);
    labels = arrayfun(@(x)sprintf('x=%.2f',x),xpos,'UniformOutput',false);
    figure; colororder(cool(length(xpos)));
    plot(real(fr),-imag(fr)); hold on;
    plot(data.dcGain,0,'o');
    plot(data.hfGain,0,'^');
    legend(labels{:},'Location','best','NumColumns',2);
    xlabel('Re');
    ylabel('-Im')
    title(varname);
    setAxesNyquist;
    thesisFormat('LineMarkerSize',markerSize);
    print(fullfile(plotdir,varname),'-depsc');
    print(fullfile(plotdir,varname),'-dpng');
end