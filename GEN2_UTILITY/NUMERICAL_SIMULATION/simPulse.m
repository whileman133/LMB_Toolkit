function simData = simPulse(varargin)
%SIMPULSE Simulate application of a charge-neutral current pulse in COMSOL.
%
% -- Usage --
% simData = simPulse(cellModel, I, Tp, socPct) simulates the
%   response of the full-order cell model CELLMODEL to a charge-neutral
%   input current pulse of magnitude I [A] (discharge followed by charge).
%   The initial SOC is set by SOCPCT [%]. TP [s] sets the duration of the
%   discharge and charge intervals (each).
%   * The COMSOL LiveLink interface for MATLAB must be running. *
%
% simData = simPulse(cellModel, I, Tp, socPct, TdegC) also
%   specifies the temperature [degC] at which the simulation is run.
%
% simData = simPulse(...,'Ts',Ts) uses sampling interval Ts [sec]
%   instead of the default (Tend/1000, where Tend is the duration of the
%   entire pulse).
%
% simData = simPulse(...,'Toffset',To) offsets the application of the pulse
%   from time t=0 by TO [sec]. This is neccesary to prevent issues with the
%   the discontinuity at t=0 in COMSOL. The default is TO=1 [sec].
%
% simData = simPulse(...,'Verbose',true) outputs status information to
%   the command window.
%
% simData = simPulse(...,'OptSimFOM',opt) passes the options structure
%   OPT to simFOM when running the COMSOL simulation.
%
% -- Output --
% The output, SIMDATA, is a structure with the following fields:
%
%   simData.time       = time vector [s]
%   simData.iapp       = applied current vector [A]
%   simData.vcell      = cell voltage vector [V]
%   simData.dis        = logical indicies to discharge interval
%   simData.chg        = logical indicies to charge interval
%   simData.tdis       = time at application of discharge pulse [s]
%   simData.tchg       = time at application of charge pulse [s] (tchg>tdis)
%   simData.Tend       = duration of the simulation [s]
%   simData.param      = structure of parameter values supplied to function
%
% -- Changelog --
% 2023.06.09 | Created | Wesley Hileman <whileman@uccs.edu>

parser = inputParser;
parser.addRequired('cellModel',@(x)isscalar(x)&&isstruct(x));
parser.addRequired('I',@(x)isscalar(x)&&x>0);
parser.addRequired('Tp',@(x)isscalar(x)&&x>0);
parser.addRequired('socPct',@(x)isscalar(x)&&0<=x&&x<=100);
parser.addOptional('TdegC',25,@(x)isscalar(x));
parser.addParameter('Ts',[],@(x)isscalar(x)&&x>0);
parser.addParameter('Toffset',1,@(x)isscalar(x)&&x>0);
parser.addParameter('Verbose',false,@(x)isscalar(x)&&islogical(x));
parser.addParameter('OptSimFOM',struct,@isstruct);
parser.parse(varargin{:});
p = parser.Results;  % structure of validated arguments

% Generate time and iapp vectors.
Tend = 2*(p.Tp+p.Toffset);
if isempty(p.Ts)
    p.Ts = Tend/1000;
end
time = 0:p.Ts:Tend;             % time vector [s]
tdis = p.Toffset;               % time at application of discharge pulse
tchg = tdis + p.Tp + p.Toffset; % time at application of charge pulse
dis = tdis<=time&time<tchg;     % logical indicies to discharge interval
chg = tchg<=time;               % logical indicies to charge interval
iapp = zeros(size(time));       % iapp vector [A]
iapp(dis) = +p.I;
iapp(chg) = -p.I;

% Generate COMSOL model.
if p.Verbose
    fprintf('Generating COMSOL model for pulse simulation... ');
end
modelCOMSOL = genFOM(p.cellModel,'DebugFlag',false);
if p.Verbose
    fprintf('done!\n');
end

% Perform simulation in COMSOL.
simspec.time = time;
simspec.Iapp = iapp;
simspec.SOC0 = p.socPct;
simspec.T = p.TdegC;
simspec.TSHIFT = 0;
if p.Verbose
    fprintf('Running pulse simulation...\n');
end
[~,sim] = simFOM(modelCOMSOL,simspec,'InputType','lut',p.OptSimFOM);
if p.Verbose
    fprintf('done!\n');
end

simData.time = time(:);
simData.iapp = iapp(:);
simData.vcell = sim.Vcell(:);
simData.chg = chg(:);
simData.dis = dis(:);
simData.tdis = tdis;
simData.tchg = tchg;
simData.Tend = Tend;
simData.param = p;
simData.origin__ = 'simPulse';

end