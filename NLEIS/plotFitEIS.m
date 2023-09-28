% plotFitEIS.m
%
% Plot regressed EIS model against lab data.

clear; close all; clc;
addpath(fullfile('..','..'));
TB.addpaths;

filename = '202309_EIS-16degC26degC-Ds=linear-k0=linear';
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
fdlp = 1./taudlp./2./pi;

% Compare lab to fit impedance: Nyquist.
fprintf('%10s%10s%10s%10s%10s\n','TdegC','MagRMSE','PhsRMSE','RealRMSE','ImagRMSE');
for k = 1:length(fitData.Zmodel)
    cellName = fitData.arg.labSpectra(k).cellName;
    cellName = split(cellName,'_');
    cellName = cellName{1};
    TdegC = TdegCvect(k);
    nsoc = length(fitData.socPctTrue{k});
    indSOC = [1 2 3:2:nsoc-2 nsoc-1 nsoc];
    nsocPlot = length(indSOC);
    freq = fitData.freq{k};
    Zlab = fitData.Zlab{k};
    Zmodel = fitData.Zmodel{k};
    rmseMag = sqrt(mean((abs(Zlab(:))-abs(Zmodel(:))).^2));
    rmsePhs = sqrt(mean((angle(Zlab(:))-angle(Zmodel(:))).^2))*180/pi;
    rmseReal = sqrt(mean((real(Zlab(:))-real(Zmodel(:))).^2));
    rmseImag = sqrt(mean((imag(Zlab(:))-imag(Zmodel(:))).^2));
    fprintf('%10.0f%10.3f%10.3f%10.3f%10.3f\n', ...
        TdegC,rmseMag,rmsePhs,rmseReal,rmseImag);

    [~,indbreakf] = min(abs(fdlp(:).'-freq(:)));
    
    figure;
    l = tiledlayout(3,ceil(nsocPlot/3));
    l.Title.String = sprintf( ...
        'Linear Impedance Spectra: Cell %s (%.0f\\circC)',cellName,TdegC);
    l.XLabel.String = 'Z'' [\Omega]';
    l.YLabel.String = '-Z'''' [\Omega]';
    for j = 1:nsocPlot
        ind = indSOC(j);
        zlab = Zlab(:,ind);
        zmod = Zmodel(:,ind);
        nexttile;
        plot(real(zlab),-imag(zlab),'b.'); hold on;
        plot(real(zmod),-imag(zmod),'r');
        plot(real(zmod(indbreakf(ind))),-imag(zmod(indbreakf(ind))),'k+');
        title(sprintf('%.0f%% SOC',fitData.socPctTrue{k}(ind)));
        setAxesNyquist;
    end
    legend('Lab','Fit Model','Location','northwest');
    thesisFormat('FigSizeInches',[10 6]);
    print(fullfile(plotdir,sprintf('Z-Nyq-%.0fdegC',TdegC)),'-depsc');
    print(fullfile(plotdir,sprintf('Z-Nyq-%.0fdegC',TdegC)),'-dpng');
end

% Plot i0.
figure;
semilogy(theta,ctData.i0,'k'); hold on;
semilogy( ...
    fitData.values.pos.k0Theta, ...
    fitData.values.pos.k0Linear,'ro');
set(gca,'xdir','reverse');
xlabel('$x$ in $\mathrm{Li}_x\mathrm{Ni}_y\mathrm{Mn}_z\mathrm{Co}_{1-y-z}\mathrm{O}_2$', ...
    'Interpreter','latex');
ylabel('Exchange current, i_{0} [A]');
title('Exchange current (linear interp.)');
thesisFormat;
print(fullfile(plotdir,'i0-theta'),'-depsc');
print(fullfile(plotdir,'i0-theta'),'-dpng');

% Plot taudlp.
figure;
semilogy(theta,taudlp,'k');
set(gca,'xdir','reverse');
xlabel('$x$ in $\mathrm{Li}_x\mathrm{Ni}_y\mathrm{Mn}_z\mathrm{Co}_{1-y-z}\mathrm{O}_2$', ...
    'Interpreter','latex');
ylabel('Double-layer time constant, R_{ct}C_{dl} [s]');
title('Double-layer time constant');
thesisFormat;
print(fullfile(plotdir,'taudlp-theta'),'-depsc');
print(fullfile(plotdir,'taudlp-theta'),'-dpng');

% Plot fbreak_dlp.
figure;
semilogy(theta,fdlp,'k');
set(gca,'xdir','reverse');
xlabel('$x$ in $\mathrm{Li}_x\mathrm{Ni}_y\mathrm{Mn}_z\mathrm{Co}_{1-y-z}\mathrm{O}_2$', ...
    'Interpreter','latex');
ylabel('Break frequency, f_b=1/2\piR_{ct}C_dl [Hz]');
title('Double-layer break frequency');
thesisFormat;
print(fullfile(plotdir,'fbreak_dlp-theta'),'-depsc');
print(fullfile(plotdir,'fbreak_dlp-theta'),'-dpng');

% Plot Ds.
figure;
semilogy(theta,dsData.Ds,'k'); hold on;
semilogy( ...
    fitData.values.pos.DsTheta, ...
    fitData.values.pos.DsLinear,'ro');
set(gca,'xdir','reverse');
xlabel('$x$ in $\mathrm{Li}_x\mathrm{Ni}_y\mathrm{Mn}_z\mathrm{Co}_{1-y-z}\mathrm{O}_2$', ...
    'Interpreter','latex');
ylabel('Solid diffusion coefficient, D_s [s^{-1}]');
title('Solid diffusion coefficient (linear interp.)');
thesisFormat;
print(fullfile(plotdir,'Ds-theta'),'-depsc');
print(fullfile(plotdir,'Ds-theta'),'-dpng');