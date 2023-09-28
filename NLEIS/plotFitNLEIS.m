% plotFitNLEIS.m
%
% Plot regressed NLEIS model against lab data.

clear; close all; clc;
addpath(fullfile('..','..'));
TB.addpaths;

filename = 'NLEIS-16degC26degC';
fitData = load(fullfile('labfitdata',[filename '.mat']));

TdegCvect = fitData.TdegC;
TrefdegC = fitData.arg.labLinearFit.arg.TrefdegC;
tmpStr = sprintf('%.0fdegC',TdegCvect);
plotdir = fullfile('plots',sprintf('LAB-NL-%s',tmpStr));
if ~isfolder(plotdir)
    mkdir(plotdir);
end

soc = linspace(0,1,100);
socPct = soc*100;
theta0 = fitData.values.pos.theta0;
theta100 = fitData.values.pos.theta100;
theta = theta0 + soc*(theta100-theta0);
electrode = MSMR(fitData.values.pos);

% Compare lab to fit impedance.
fprintf('%10s%10s%10s%10s%10s\n','TdegC','MagRMSE','PhsRMSE','RealRMSE','ImagRMSE');
for k = 1:length(fitData.Zmodel)
    cellName = fitData.arg.labLinearFit.arg.labSpectra(k).cellName;
    cellName = split(cellName,'_');
    cellName = cellName{1};
    TdegC = TdegCvect(k);
    nsoc = length(fitData.socPctTrue{k});
    indSOC = [1 2 3:2:nsoc-2 nsoc-1 nsoc];
    nsocPlot = length(indSOC);
    Zlab = fitData.Zlab{k};
    Zmodel = fitData.Zmodel{k};
    freq = fitData.freq{k};
    rmseMag = sqrt(mean((abs(Zlab(:))-abs(Zmodel(:))).^2));
    rmsePhs = sqrt(mean((angle(Zlab(:))-angle(Zmodel(:))).^2))*180/pi;
    rmseReal = sqrt(mean((real(Zlab(:))-real(Zmodel(:))).^2));
    rmseImag = sqrt(mean((imag(Zlab(:))-imag(Zmodel(:))).^2));
    fprintf('%10.0f%10.3f%10.3f%10.3f%10.3f\n', ...
        TdegC,rmseMag,rmsePhs,rmseReal,rmseImag);
    
    % Nyquist.
    figure;
    l = tiledlayout(3,ceil(nsocPlot/3));
    l.Title.String = sprintf( ...
        'Second-Harmonic Impedance Spectra: Cell %s (%.0f\\circC)',cellName,TdegC);
    l.XLabel.String = 'Z'' [\Omega A^{-1}]';
    l.YLabel.String = '-Z'''' [\Omega A^{-1}]';
    for j = 1:nsocPlot
        ind = indSOC(j);
        zlab = Zlab(:,ind);
        zmod = Zmodel(:,ind);
        nexttile;
        plot(real(zmod),-imag(zmod),'r');  hold on;
        plot(real(zlab),-imag(zlab),'b.');
        title(sprintf('%.0f%% SOC',fitData.socPctTrue{k}(ind)));
        setAxesNyquist;
    end
    legend('Model','Lab','Location','best');
    thesisFormat('FigSizeInches',[10 6]);
    print(fullfile(plotdir,sprintf('Z-Nyq-%.0fdegC',TdegC)),'-depsc');
    print(fullfile(plotdir,sprintf('Z-Nyq-%.0fdegC',TdegC)),'-dpng');

    % Bode Magnitude.
    figure;
    l = tiledlayout(3,ceil(nsocPlot/3));
    l.Title.String = sprintf( ...
        'Second-Harmonic Impedance Spectra: Cell %s (%.0f\\circC)',cellName,TdegC);
    l.XLabel.String = 'Cyclic Frequency, f [Hz]';
    l.YLabel.String = '|Z| [\Omega A^{-1}]';
    for j = 1:nsocPlot
        ind = indSOC(j);
        zlab = Zlab(:,ind);
        zmod = Zmodel(:,ind);
        nexttile;
        loglog(freq,abs(zmod),'r');  hold on;
        loglog(freq,abs(zlab),'b.');
        xlim([min(freq) max(freq)]);
        title(sprintf('%.0f%% SOC',fitData.socPctTrue{k}(ind)));
    end
    legend('Model','Lab','Location','best');
    thesisFormat('FigSizeInches',[10 6]);
    print(fullfile(plotdir,sprintf('Z-BodeMag-%.0fdegC',TdegC)),'-depsc');
    print(fullfile(plotdir,sprintf('Z-BodeMag-%.0fdegC',TdegC)),'-dpng');

    % Bode Phase.
    figure;
    l = tiledlayout(3,ceil(nsocPlot/3));
    l.Title.String = sprintf( ...
        'Second-Harmonic Impedance Spectra: Cell %s (%.0f\\circC)',cellName,TdegC);
    l.XLabel.String = 'Cyclic Frequency, f [Hz]';
    l.YLabel.String = '\angleZ [deg]';
    for j = 1:nsocPlot
        ind = indSOC(j);
        zlab = Zlab(:,ind);
        zmod = Zmodel(:,ind);
        nexttile;
        semilogx(freq,angle(zmod)*180/pi,'r');  hold on;
        semilogx(freq,angle(zlab)*180/pi,'b.');
        xlim([min(freq) max(freq)]);
        title(sprintf('%.0f%% SOC',fitData.socPctTrue{k}(ind)));
    end
    legend('Model','Lab','Location','best');
    thesisFormat('FigSizeInches',[10 6]);
    print(fullfile(plotdir,sprintf('Z-BodePhs-%.0fdegC',TdegC)),'-depsc');
    print(fullfile(plotdir,sprintf('Z-BodePhs-%.0fdegC',TdegC)),'-dpng');
end