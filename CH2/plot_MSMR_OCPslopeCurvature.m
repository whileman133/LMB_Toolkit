clear; close all; clc;

if ~exist('assets','dir')
    % Bootstrap the toolbox.
    addpath('..');
    TB.addPaths();
end

el = MSMR.C6();
figs = el.plotOCP();

exportgraphics(figs.ocp,'MSMR-ocp.eps');
exportgraphics(figs.docp,'MSMR-docp.eps');
exportgraphics(figs.d2ocp,'MSMR-d2ocp.eps');
