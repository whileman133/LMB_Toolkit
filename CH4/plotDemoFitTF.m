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
