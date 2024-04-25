clear; close all; clc;
addpath('..');
TB.addpaths('gen2');

el = MSMR.NMC622();
figs = el.plotOCP('vmin',3.2,'vmax',4.3,'npoints',100);
figure(figs.ocp);
ylim([3.2 4.3]);

exportgraphics(figs.ocp,'NMC622-ocp.eps');
exportgraphics(figs.ocp,'NMC622-ocp.png');
exportgraphics(figs.docp,'NMC622-docp.eps');
exportgraphics(figs.d2ocp,'NMC622-d2ocp.eps');
