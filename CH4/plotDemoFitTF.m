% plotDemoFitTF.m

clear; close all; clc;
addpath('..');
TB.addpaths;
load('demoFitTF.mat');
fitmodel = fitdata.values;
modelspec = fitdata.modelspec;

% Evalulate parameter estimates at reference temperature.
TrefdegC = getCellParams(rlwrm,'const.Tref')-273.15;
sparse = fastopt.unpack(fastopt.pack(fitmodel,modelspec),modelspec,'sparse',true);
estimates = fastopt.evaltemp(sparse,modelspec,TrefdegC);
estimates = orderfields(estimates,{'const','neg','eff','pos','pkg'});

% Evalulate truth at reference temperature.
truth = getCellParams(rlwrm,'TdegC',TrefdegC);
truth.pkg.R0 = truth.pkg.R0 + truth.neg.Rf;  % fit model lumps Rfn with R0

% Compare estimates to truth.
clear names est tru pctErr lb ub EactKJ;
cursor = 1;
rnames = fieldnames(estimates);
for kr = 1:length(rnames)
    rname = rnames{kr};
    pnames = fieldnames(estimates.(rname));
    for kp = 1:length(pnames)
        pname = pnames{kp};
        e = estimates.(rname).(pname);
        t = truth.(rname).(pname);
        lwr = fitdata.lb.(rname).(pname);
        upr = fitdata.ub.(rname).(pname);
        m = length(e);
        for km = 1:m
            if m == 1
                names(cursor) = sprintf("%s.%s",rname,pname);
            else
                names(cursor) = sprintf("%s.%s_%d",rname,pname,km);
            end
            est(cursor) = e(km);
            tru(cursor) = t(km);
            pctErr(cursor) = 100*(est(cursor)-tru(cursor))./tru(cursor);
            lb(cursor) = lwr(km);
            ub(cursor) = upr(km);
            if isfield(sparse.(rname),[pname '_Eact'])
                EactKJ(cursor) = sparse.(rname).([pname '_Eact'])/1000;
            else
                EactKJ(cursor) = NaN;
            end
            cursor = cursor + 1;
        end
    end
end
names = names.';
est = est.';
tru = tru.';
pctErr = pctErr.';
lb = lb.';
ub = ub.';
EactKJ = EactKJ.';
tab = table;
tab.Name = names;
tab.Estimate = est;
tab.Truth = tru;
tab.PctError = pctErr;
tab.LowerBound = lb;
tab.UpperBound = ub;
tab.EactKJ = EactKJ;

% Save results to spreadsheet.
writetable(tab,'demoFitTF.xlsx');

% Plotting ----------------------------------------------------------------
TdegC = 25;
lab = fitdata.lab;
fit = fitdata.fit;
indT = [lab.TdegC]==TdegC;
lab = lab(indT);
fit = fit(indT);
nsoc = length(lab.socPct);
labels = arrayfun( ...
    @(kz)sprintf('%.0f%%',lab.socPct(kz)),1:nsoc, ...
    'UniformOutput',false);

figure; colororder(cool(nsoc));
plot(real(fit.Z),-imag(fit.Z)); hold on;
plot(real(lab.Z),-imag(lab.Z),'o');
xlabel('Z_{cell}'' [\Omega]');
ylabel('-Z_{cell}'''' [\Omega]');
title('Nyquist: Z_{cell} (T=25\circC)');
legend(labels{:},'NumColumns',4,'Location','northwest','FontSize',11);
setAxesNyquist('ydata',[0 0.5]);
thesisFormat('LineMarkerSize',2,'LineMarkerLineWidth',0.5,'LineLineWidth',1);
ax = addInset([3.5 3.65],[3.25 0.12],2.3,'YSpan',[0 0.08]);
setAxesNyquist('axes',ax,'xdata',[3.5 3.65],'ydata',[0 0.08]);
print(fullfile('plots',sprintf('Z-Nyq-SINGLE-%.0fdegC',TdegC)),'-depsc');
print(fullfile('plots',sprintf('Z-Nyq-SINGLE-%.0fdegC',TdegC)),'-dpng');

