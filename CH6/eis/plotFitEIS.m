% plotFitEIS.m
%
% Plot regressed EIS model against lab data.

clear; close all; clc;
addpath(fullfile('..','..'));
TB.addpaths;

filename = 'EIS-16degC26degC-Ds=linear-k0=linear';
fitData = load(fullfile('labfitdata',[filename '.mat']));

TdegCvect = fitData.TdegC;
TrefdegC = fitData.arg.TrefdegC;
tmpStr = sprintf('%.0fdegC',TdegCvect);
plotdir = fullfile('plots',sprintf('LAB-%s',tmpStr));
if ~isfolder(plotdir)
    mkdir(plotdir);
end

soc = linspace(0,1,100);
socPct = soc*100;
theta0 = fitData.values.pos.theta0;
theta100 = fitData.values.pos.theta100;
theta = theta0 + soc*(theta100-theta0);
electrode = MSMR(fitData.values.pos);
ctData = electrode.Rct(fitData.values.pos,'TdegC',TrefdegC,'theta',theta);
dsData = electrode.Ds(fitData.values.pos,'TdegC',TrefdegC,'theta',theta);
Cdlp = fitData.values.pos.Cdl;
taudlp = Cdlp.*ctData.Rct;

% Compare lab to fit impedance: Nyquist.
fprintf('%10s%10s%10s%10s%10s\n','TdegC','MagRMSE','PhsRMSE','RealRMSE','ImagRMSE');
for k = 1:length(fitData.Zmodel)
    cellName = fitData.arg.labSpectra(k).cellName;
    cellName = split(cellName,'_');
    cellName = cellName{1};
    TdegC = TdegCvect(k);
    nsoc = length(fitData.socPctTrue{k});
    freq = fitData.freq{k};
    Zlab = fitData.Zlab{k};
    Zmodel = fitData.Zmodel{k};
    rmseMag = sqrt(mean((abs(Zlab(:))-abs(Zmodel(:))).^2));
    rmsePhs = sqrt(mean((angle(Zlab(:))-angle(Zmodel(:))).^2))*180/pi;
    rmseReal = sqrt(mean((real(Zlab(:))-real(Zmodel(:))).^2));
    rmseImag = sqrt(mean((imag(Zlab(:))-imag(Zmodel(:))).^2));
    fprintf('%10.0f%10.3f%10.3f%10.3f%10.3f\n', ...
        TdegC,rmseMag*1000,rmsePhs,rmseReal*1000,rmseImag*1000);

    lab = arrayfun( ...
        @(kz)sprintf('%.0f%%',fitData.socPctTrue{k}(kz)),1:nsoc, ...
        'UniformOutput',false);
    figure; colororder(cool(nsoc));
    plot(real(Zmodel),-imag(Zmodel)); hold on;
    plot(real(Zlab),-imag(Zlab),'o');
    xlabel('re');
    ylabel('-im');
    title('nyquist');
    legend(lab{:},'NumColumns',4,'Location','southeast','FontSize',8);
    thesisFormat('LineMarkerSize',2,'LineMarkerLineWidth',0.5,'LineLineWidth',1);
    if abs(TdegC-25)<=2
        ax = addInset([0.1 1.4],[3.5 4.5],2,'YSpan',[0 0.8]);
        setAxesNyquist('axes',ax,'xdata',[0.1 1.4],'ydata',[0 0.8]);
    else
        ax = addInset([0.1 1.4],[5 4.5],2,'YSpan',[0 0.8]);
        setAxesNyquist('axes',ax,'xdata',[0.1 1.4],'ydata',[0 0.8]);
    end
    print(fullfile(plotdir,sprintf('Z-Nyq-SINGLE-%.0fdegC',TdegC)),'-depsc');
    print(fullfile(plotdir,sprintf('Z-Nyq-SINGLE-%.0fdegC',TdegC)),'-dpng');
end