% plotOCPRaw.m
%
% Plot OCP data collected during the four scripts and the calibrated
% dis/charge curves.
%
% -- Changelog --
% 10.12.2022 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath('..');

% Constants.
studyname = 'SionFresh_0C01';
TdegC = 40;
plotdir = fullfile('plots',sprintf('RAW_%s_%.0fdegC',studyname,TdegC));

ocp.load(studyname,'built')
test = builtstudy.getTest(TdegC);

if ~isfolder(plotdir)
    mkdir(plotdir);
end

figure;
plot(test.script1.time/3600,test.script1.voltage);
xlabel('time');
ylabel('vcell');
title('Script1');
thesisFormat;
print(fullfile(plotdir,'script1'),'-depsc');
print(fullfile(plotdir,'script1'),'-dpng');

figure;
plot(test.script2.time/3600,test.script2.voltage);
xlabel('time');
ylabel('vcell');
title('Script2');
thesisFormat;
addInset([0.91 1.4],[0.7 3.3]);
print(fullfile(plotdir,'script2'),'-depsc');
print(fullfile(plotdir,'script2'),'-dpng');

figure;
plot(test.script3.time/3600,test.script3.voltage);
xlabel('time');
ylabel('vcell');
title('Script3');
thesisFormat;
print(fullfile(plotdir,'script3'),'-depsc');
print(fullfile(plotdir,'script3'),'-dpng');

figure;
plot(test.script4.time/3600,test.script4.voltage);
ylim([4.28 4.32]);
xlabel('time');
ylabel('vcell');
title('Script4');
thesisFormat;
addInset([0.91 1.4],[0.7 4.304]);
print(fullfile(plotdir,'script4'),'-depsc');
print(fullfile(plotdir,'script4'),'-dpng');

figure;
plot(test.disZ,test.disV,'-'); hold on;
plot(test.chgZ,test.chgV,':');
xlabel('SOC');
ylabel('vcell');
title('CalSOC');
legend('Discharge','Charge');
thesisFormat;
print(fullfile(plotdir,'VdisVchg-v-SOC'),'-depsc');
print(fullfile(plotdir,'VdisVchg-v-SOC'),'-dpng');