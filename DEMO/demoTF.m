% demoTF.m
%
% Demonstrate the usage of the toolkit for computing transfer functions 
% of electrochemical variables for lithium-metal battery cells.
%
% NOTE: All of the TF implementation appears in UTILITY/CELL_MODEL/tfLMB.m
% The tfXX.m functions call the tfLMB.m function internally.
%
% NOTE: For half cells, use the "effective" layer as the separator 
% (eff => sep).
%
% -- Changelog --
% 2024.01.18 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath('..');
TB.addpaths;

% Constants.
cellsheet = 'cellNMC30-RLP2DM.xlsx';  % spreadsheet with param values
freq = logspace(-3.5,5,100);          % vector of frequency points [Hz]
socPct = 10;  % soc at which to evalulate impedance [%]
TdegC  = 25;  % temperature at which to evalulate impedance [degC]

% Load spreadsheet of cell parameter values.
cellmodel = loadCellModel(cellsheet);

% Convert cell model to 'Legacy Lumped-Parameter Model' form, which
% works with Plett's toolbox for lithium-ion batteries.
cellmodel = convertCellModel( ...
    cellmodel,'LLPM', ...
    'LegacyExpandEff',false ... % !!! important, otherwise x goes all the way to 3
);

% Compute impedance of the cell using tfXX functions:
s = 1j*2*pi*freq;
cellparams = evalSetpoint(cellmodel,s,socPct/100,TdegC+273.15);
Phise  = tfPhiseInt(s,[0 2],cellparams);
Phise0 = Phise(1,:);  % Phise/Iapp impedance at x=0
Phise2 = Phise(2,:);  % Phise/Iapp impedance at x=2 (pos current collector)
Phie2  = tfPhie(s,2,cellparams);  % Phie/Iapp impedance at x=2 (pos current collector)
Vcell  = -Phise0 + Phie2 + Phise2; % Vcell/Iapp
Zcell  = -Vcell; % -Vcell/Iapp is the cell impedance

% Or, equivalently, use tfLMB.m directly:
%tf = tfLMB(s,cellmodel,'TdegC',TdegC,'socPct',socPct);
%Zcell = tf.h11.tfZcell();  % h11 indicates the linear (first harmonic)

% Plot (Nyquist).
figure;
plot(real(Zcell),-imag(Zcell));
xlabel("Z' [\Omega]");
ylabel("-Z'' [\Omega]");
title(cellmodel.name);
if exist('quadprog','file')
    % Need quadprog to format the axes limits so that the Nyquist
    % plot is true 1:1 scale
    setAxesNyquist;
end
thesisFormat;