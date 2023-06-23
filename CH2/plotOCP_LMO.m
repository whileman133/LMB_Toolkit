clear; close all; clc;
addpath('..');
TB.addpaths('gen2');

el = MSMR.LMO();
figs = el.plotOCP();

exportgraphics(figs.ocp,'LMO-ocp.eps');
exportgraphics(figs.docp,'LMO-docp.eps');
exportgraphics(figs.d2ocp,'LMO-d2ocp.eps');
