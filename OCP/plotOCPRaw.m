% plotOCPRaw.m
%
% Plot OCP data collected during the four scripts and the calibrated
% dis/charge curves.
%
% -- Changelog --
% 10.12.2022 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath('..');
ocp.load('SionFresh_0C01','built');

TdegC = 40;
test = builtstudy.getTest(TdegC);

figure;
plot(test.script1.time/3600,test.script1.voltage);
xlabel('Time [hours]');
ylabel('Terminal Voltage [V]');
title( ...
    sprintf('Script 1: C/%d Discharge (T=%dDegC)', ...
    1/builtstudy.crate,test.temp));
thesisFormat;

figure;
plot(test.script2.time/3600,test.script2.voltage);
xlabel('Time [hours]');
ylabel('Terminal Voltage [V]');
title( ...
    sprintf('Script 2: Dither at Vmin (T=%dDegC)', ...
    test.temp));
thesisFormat;

figure;
plot(test.script3.time/3600,test.script3.voltage);
xlabel('Time [hours]');
ylabel('Terminal Voltage [V]');
title( ...
    sprintf('Script 3: C/%d Charge (T=%dDegC)', ...
    1/builtstudy.crate,test.temp));
thesisFormat;

figure;
plot(test.script4.time/3600,test.script4.voltage);
xlabel('Time [hours]');
ylabel('Terminal Voltage [V]');
title( ...
    sprintf('Script 4: Dither at Vmax (T=%dDegC)', ...
    test.temp));
thesisFormat;

figure;
plot(test.disZ,test.disV); hold on;
plot(test.chgZ,test.chgV);
xlabel('Cell State-of-Charge [%]');
ylabel('Terminal Voltage');
legend('Discharge','Charge');
title('Calibrated Dis/charge Curves');
thesisFormat;