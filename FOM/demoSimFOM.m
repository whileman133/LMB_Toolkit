% demoSimFOM.m
%
% Demonstate usage of the genFOM and simFOM utility functions for
% programmatically generating and running full-order COMSOL simulations.
% Run this file in COMSOL with MATLAB LiveLink.
%
% For further information, use:
%    help genFOM
%    help simFOM
% in the command window.
%
% 2023.07.24 | Update for gen2 toolkit | Wesley Hileman
% 2023.04.10 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath('..');
TB.addpaths;

% Constants.
cellFile = 'cellSionGuess-P2DM.xlsx';  % Name of cell parameters spreadsheet.


% 1. Working with cell models. --------------------------------------------

% Load cell model from spreadsheet.
cellModel = loadCellModel(cellFile);
% Convert standard cell model to lumped-parameter model.
cellModel = convertCellModel(cellModel,'LLPM');

% `cellModel.function` is now a structure containing the parameters of the cell
% as functions of lithiation (x) and absolute temperature (T). The
% structure has fields for each region of the cell (neg,dll,sep,pos) as well
% as constants (const).

% Evaluate some cell parameters at 25C...
Q = cellModel.function.const.Q([],25+273.15);  % ok that x empty since we know Q invariant with lithiation
k0n = cellModel.function.neg.k0([],25+273.15);
k0p = cellModel.function.pos.k0([],25+273.15);

% Or evaluate all cell parameters at a particular SOC setpoint and
% temperature using evalSetpoint...
socPct = 50;
TdegC = 25;
cellModelSetpoint1 = evalSetpoint(cellModel,[],socPct/100,TdegC+273.15);

% Now, cellModelSetpoint1.neg, ".dll, ".sep, ".pos, ".const are structures
% containing the parameters evalulated at socPct and TdegC.


% 2. Generating COMSOL models. --------------------------------------------

% Create COMSOL model from cell model. LiveLink must be running!
% genData.FOM is the COMSOL model object.
genData = genFOM(cellModel);

% If you don't want the status messages output to the console, use:
% genData = genFOM(cellModel,'DebugFlag',false);

% Save standard COMSOL metaphysics (mph) file that you can open
% in the COMSOL GUI:
genData.FOM.save('demoFOM.mph');


% 3. Running a simulation in COMSOL. --------------------------------------

% First, generate input current iapp(t) and temperature T(t) waveforms. 
% You can use loadInput.m to load a profile from a spreadsheet in
% XLSX_INPUTS, or create an input programmatically in the code.

% Using loadInput...
inputUDDS = loadInput('inputUDDS.xlsx');

% Programatically...
inputCCDis.SOC0 = 100; % starting SOC [%].
inputCCDis.time = 0:0.1:100; % time vector [s].
% Applied current vector [A]. 1C constant-current dischange for this example.
inputCCDis.Iapp = zeros(size(inputCCDis.time));
inputCCDis.Iapp(inputCCDis.time>=1) = 1*Q;
inputCCDis.T = 25; % Temperature in degC. (Can also be a vector.)
inputCCDis.TSHIFT = 0; % see comment in simFOM.m, need nonzero value if Iapp(0)~=0!

% (The default input type for simFOM is 'interpolation' [aka lookup-table]
% and that is what loadInput.m supports and what we demonstrate here.
% Also available are 'sine' and 'PWM' input types. Type `help simFOM` for
% more information on using these types of inputs.)

% Run simulation in COMSOL. (This may take several seconds/minutes.)
[FOMwithResults,simData] = simFOM(genData,inputCCDis);

% simData now contains the time-domain results of the simulation. 
% We plot some quantities below...

% Vcell(t)
figure;
plot(simData.time,simData.Vcell);
xlabel('t [s]');
ylabel('V_{cell} [V]');
title('Cell voltage vs. time');
thesisFormat;

% Thetae(x,t)
% simData.Thetae is a matrix:
%  - dim1 (rows) = time
%  - dim2 (columns) = unitless x-location in cell
xxThetae = simData.xLocs.Thetae(:)';
xxDesiredPlot = 0:3;
labels = arrayfun(@(x)sprintf('x=%d',x),xxDesiredPlot,'UniformOutput',false);
[~,indxx] = min(abs(xxDesiredPlot-xxThetae'));
figure;
plot(simData.time,simData.Thetae(:,indxx));
ylabel('\theta_e(x,t) [unitless]');
xlabel('t [s]');
title('Electrolyte concentration vs time');
legend(labels{:},'Location','best');
thesisFormat;