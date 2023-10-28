function simData = simCC(varargin)
%SIMCC Simulate constant-current dis/charge of a LMB cell in COMSOL.
%
% -- Usage --
% simData = simCC(cellModel, I, soc0Pct, socfPct) simulates the
%   response of the full-order cell model CELLMODEL to a constant 
%   dis/charge current of value I [A]. The initial SOC is set by 
%   SOC0PCT [%] and the final SOC by SOCFPCT [%].
%   * The COMSOL LiveLink interface for MATLAB must be running. *
%
% simData = simCC(cellModel, Iavg, soc0Pct, socfPct, TdegC) also
%   specifies the temperature [degC] at which the simulation is run.
%
% simData = simCC(...,'Ts',Ts) uses sampling interval Ts [sec]
%   instead of the default (Tend/1000, where Tend is the duration of the
%   input).
%
% simData = simCC(...,'Trelax',Tr) allows the cell to relax for Tr [sec]
%   after the current is removed before ending the simulation.
%
% simData = simCC(...,'Verbose',true) outputs status information to
%   the command window.
%
% simData = simCC(...,'OptSimFOM',opt) passes the options structure
%   OPT to simFOM when running the COMSOL simulation.
%
% simData = simCC(...,'DryRun',true) performs all calculations, but
%   does not run the simulation in COMSOL. Useful for generating the
%   current waveform for the constant-current discharge.
%
% -- Output --
% The output, SIMDATA, is a structure with the following fields:
%
%   simData.time    = time vector [s]
%   simData.iapp    = applied current vector [A]
%   simData.vcell   = cell voltage vector [V] (omitted for 'DryRun')
%   simData.output  = data returned by simFOM (omitted for 'DryRun')
%   simData.Qdis    = total charge removed from cell (- for added) [Ah]
%   simData.Tend    = duration of the discharge waveform [s]
%   simData.arg     = structure of parameter values supplied to function
%
% -- Changelog --
% 2023.06.23 | Created | Wesley Hileman <whileman@uccs.edu>

parser = inputParser;
parser.addRequired('cellModel',@(x)isscalar(x)&&isstruct(x));
parser.addRequired('I',@(x)isscalar(x)&&x>0);
parser.addRequired('soc0Pct',@(x)isscalar(x)&&0<=x&&x<=100);
parser.addRequired('socfPct',@(x)isscalar(x)&&0<=x&&x<=100);
parser.addOptional('TdegC',25,@(x)isscalar(x));
parser.addParameter('Ts',[],@(x)isscalar(x)&&x>0);
parser.addParameter('Trelax',0,@(x)isscalar(x)&&x>=0);
parser.addParameter('Verbose',false,@(x)isscalar(x)&&islogical(x));
parser.addParameter('OptSimFOM',struct,@isstruct);
parser.addParameter('DryRun',false,@(x)isscalar(x)&&islogical(x));
parser.parse(varargin{:});
arg = parser.Results;  % structure of validated params

% Compute cell capacity to discharge (will be negative for charge).
Qtot = getCellParams(arg.cellModel,'const.Q');
Qdis = Qtot*(arg.soc0Pct-arg.socfPct)/100;   % capacity to discharge [Ah] 

% Compute constant-current waveform parameters.
Tend = 3600*abs(Qdis)/arg.I;    % time at end of dis/charge [s]
Tsim = Tend+arg.Trelax;         % time at end of simulation [s]
if isempty(arg.Ts)
    arg.Ts = Tend/1000;         % default sampling interval [s]
end
time = 0:arg.Ts:Tsim;           % time vector [s]
iapp = sign(Qdis)*abs(arg.I)*ones(size(time));  % applied current vector [A]
iapp(time>Tend) = 0;

if ~arg.DryRun
    % Generate COMSOL model.
    if arg.Verbose
        fprintf('Generating COMSOL model for CC simulation... ');
    end
    modelCOMSOL = genFOM(arg.cellModel,'DebugFlag',false);
    if arg.Verbose
        fprintf('done!\n');
    end

    % Perform simulation in COMSOL.
    simspec.time = time;
    simspec.Iapp = iapp;
    simspec.SOC0 = arg.soc0Pct;
    simspec.T = arg.TdegC;
    simspec.TSHIFT = arg.Ts;
    if arg.Verbose
        fprintf('Running CC simulation...\n');
    end
    [~,sim] = simFOM(modelCOMSOL,simspec,'InputType','lut',arg.OptSimFOM);
    if arg.Verbose
        fprintf('done!\n');
    end
end

simData.time = time(:);
simData.iapp = iapp(:);
if ~arg.DryRun
    simData.vcell = sim.Vcell(:);
    simData.output = sim;
end
simData.Qdis = Qdis;
simData.Tend = Tend;
simData.Tsim = Tsim;
simData.arg = arg;
simData.origin__ = 'simCC';

end