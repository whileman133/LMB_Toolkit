% plotOverlithiation.m

clear; close all; clc;
addpath('..');
TB.addpaths;

data = loadGamryDTA(fullfile('labdata','OverlithiationDischarge_1V.DTA'));
iappDis = -data.tables.CURVE.Im;
vcellDis = data.tables.CURVE.Vf;
timeDis = data.tables.CURVE.T;
QdisAh = cumtrapz(timeDis,iappDis)/3600;
figure;
plot(QdisAh,vcellDis);
xlabel('Charge Removed [Ah]');
ylabel('Cell Voltage [V]');
title('Overlithiation Discharge to ~1V [iapp=5.8mA]');
thesisFormat;
print(fullfile('plots','overlithiation-discharge-1v'),'-dpng');

data = loadGamryDTA(fullfile('labdata','OverlithiationDischarge_0V.DTA'));
iappDis = -data.tables.CURVE.Im;
vcellDis = data.tables.CURVE.Vf;
timeDis = data.tables.CURVE.T;
QdisAh = cumtrapz(timeDis,iappDis)/3600;
figure;
plot(QdisAh,vcellDis);
xlabel('Charge Removed [Ah]');
ylabel('Cell Voltage [V]');
title('Overlithiation Discharge to 0V [iapp=5.8mA]');
thesisFormat;
print(fullfile('plots','overlithiation-discharge-0v'),'-dpng');

