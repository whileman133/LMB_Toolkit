clear; close all; clc;
addpath('..');
TB.addpaths('gen2');

el = MSMR.C6();
figs = el.plotOCP();

exportgraphics(figs.ocp,'MSMR-ocp.eps');
exportgraphics(figs.docp,'MSMR-docp.eps');
exportgraphics(figs.d2ocp,'MSMR-d2ocp.eps');
