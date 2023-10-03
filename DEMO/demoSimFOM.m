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
% There are also several functions available to simulate specific types of
% inputs. Use the help command on any of the following for more info.
%
%   simCC.m: simulate constant-current discharge
%   simPulse.m: simulate charge-neutral pulse
%   sumHalfCycle.m: simulate half-cycle discharge
%   simEIS.m: simulate EIS at given SOC and frequency points. Process the
%     collected data into impedance spectra using processEIS.m
%
% 2023.07.24 | Update for gen2 toolkit | Wesley Hileman
% 2023.04.10 | Created | Wesley Hileman <whileman@uccs.edu>

clear; close all; clc;
addpath('..');
TB.addpaths;

% Constants.
cellFile = 'cellLMO-P2DM.xlsx';  % Name of cell parameters spreadsheet.


% 1. Working with cell models. --------------------------------------------

% Load cell model from spreadsheet.
p2dm = loadCellModel(cellFile);  % pseudo two-dimensional model
% Convert standard cell model to lumped-parameter Warburg-resistance model.
wrm = convertCellModel(p2dm,'WRM');

% Evaluate some cell parameters at 25degC...
% Use `help getCellParams.m` for more information.
[Q,k0n,k0p] = getCellParams(wrm,'const.Q *.k0','TdegC',25,'Output','list');

% Or evaluate all cell parameters at a particular SOC setpoint and
% temperature...
cellModelSetpoint1 = getCellParams(wrm,'TdegC',25,'socPct',50);



% 2. Generating COMSOL models. --------------------------------------------

% Create COMSOL model from cell model. LiveLink must be running!
% genData.FOM is the COMSOL model object.
genData = genFOM(wrm);

% If you don't want the status messages output to the console, use:
% genData = genFOM(wrm,'DebugFlag',false);

% Save standard COMSOL metaphysics (mph) file that you can open
% in the COMSOL GUI:
genData.FOM.save('demoFOM.mph');


% 3. Running a simulation in COMSOL. --------------------------------------

% First, generate input current iapp(t) and temperature T(t) waveforms. 
% You can use loadInput.m to load a profile from a spreadsheet in
% XLSX_INPUTS, or create an input programmatically in the code.

% Using loadInput...
inputUDDS = loadInput('inputUDDS.xlsx');
% Rescale avg dis/charge current to 1C.
Iavg = mean(abs(inputUDDS.Iapp));
inputUDDS.Iapp = inputUDDS.Iapp*(Q/Iavg);
% Ensure SOC bounds are not violated.
QdisAh = cumtrapz(inputUDDS.time,inputUDDS.Iapp);
socAvgPct = inputUDDS.SOC0 - 100*QdisAh/3600/Q;
if max(socAvgPct)>95 || min(socAvgPct)<5
    error('UDDS cycle brings cell too close to SOC limits.')
end

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
[FOMwithResults,simData] = simFOM(genData,inputUDDS);

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