function simData = simHalfCycle(varargin)
%SIMHALFCYCLE Simulate half-cycle dis/charge of a LMB cell in COMSOL. Useful
%  for comparison to the perturbation-resistance approximation developed by 
%  Baker and Verbrugge [1].
%
% -- Usage --
% simData = simHalfCycle(cellModel, Iavg, soc0Pct, socfPct) simulates the
%   response of the full-order cell model CELLMODEL to a sinusoidal
%   half-wave dis/charge current with average value IAVG [A]. The initial 
%   SOC is set by SOC0PCT [%] and the final SOC by SOCFPCT [%]. 
%   * The COMSOL LiveLink interface for MATLAB must be running. *
%
% simData = simHalfCycle(cellModel, Iavg, soc0Pct, socfPct, TdegC) also
%   specifies the temperature [degC] at which the simulation is run.
%
% simData = simHalfCycle(...,'Ts',Ts) uses sampling interval Ts [sec]
%   instead of the default (Tend/1000, where Tend is the duration of the
%   half-wave).
%
% simData = simHalfCycle(...,'Verbose',true) outputs status information to
%   the command window.
%
% simData = simHalfCycle(...,'OptSimFOM',opt) passes the options structure
%   OPT to simFOM when running the COMSOL simulation.
%
% simData = simHalfCycle(...,'DryRun',true) performs all calculations, but
%   does not run the simulation in COMSOL. Useful for generating the
%   current waveform for the half-cycle discharge
%
% -- Output --
% The output, SIMDATA, is a structure with the following fields:
%
%   simData.time       = time vector [s]
%   simData.iapp       = applied current vector [A]
%   simData.vcell      = cell voltage vector [V] (omitted for 'DryRun')
%   simData.zavg       = average fractional SOC vector (porous electrode)
%   simData.thetaAvg   = average stiochiometry vector (porous electrode)
%   simData.Qdis       = total charge removed from cell (- for added) [Ah]
%   simData.Ipk        = peak value (amplitude) of the iapp half-wave [A]
%   simData.f0         = cyclic frequency of the iapp sinusoid [Hz]
%   simData.Tend       = duration of the half-wave [s]
%   simData.param      = structure of parameter values supplied to function
%
% -- References --
% [1] Daniel R. Baker and Mark W. Verbrugge 2021 J. Electrochem. Soc. 168 050526
%
% -- Changelog --
% 2023.06.06 | Add 'DryRun' option | Wesley Hileman
% 2023.05.30 | Created | Wesley Hileman <whileman@uccs.edu>

parser = inputParser;
parser.addRequired('cellModelOrQ',@(x)isscalar(x)&&(isstruct(x)||isnumeric(x)));
parser.addRequired('Iavg',@(x)isscalar(x)&&x>0);
parser.addRequired('soc0Pct',@(x)isscalar(x)&&0<=x&&x<=100);
parser.addRequired('socfPct',@(x)isscalar(x)&&0<=x&&x<=100);
parser.addOptional('TdegC',25,@(x)isscalar(x));
parser.addParameter('Ts',[],@(x)isscalar(x)&&x>0);
parser.addParameter('Verbose',false,@(x)isscalar(x)&&islogical(x));
parser.addParameter('OptSimFOM',struct,@isstruct);
parser.addParameter('DryRun',false,@(x)isscalar(x)&&islogical(x));
parser.parse(varargin{:});
p = parser.Results;  % structure of validated params

cellModelType = getCellModelType(p.cellModelOrQ,false);
if isempty(cellModelType)
    Q = p.cellModelOrQ;  % not a cell model, interpret as capacity!
    theta0 = NaN;
    theta100 = NaN;
else
    cellModel = convertCellModel(p.cellModelOrQ,'LLPM');
    [Q,theta0,theta100] = getCellParams(cellModel, ...
        'const.Q pos.theta0 pos.theta100','Output','list');
end

% Compute cell capacity to discharge (will be negative for charge).
Qdis = Q*(p.soc0Pct-p.socfPct)/100;   % capacity to discharge [Ah] 

% Compute half-cycle waveform parameters.
Ipk = sign(Qdis)*p.Iavg*pi/2; % amplitude of iapp half-sinusoid [A]
f0 = p.Iavg/abs(Qdis)/3600/2; % cyclic frequency of iapp sinusoid [Hz]
Tend = 1/f0/2;                % time at end of half-cycle [s]
if isempty(p.Ts)
    p.Ts = Tend/1000;         % default sampling interval [s]
end
time = 0:p.Ts:Tend;           % time vector [s]
zavg = (p.soc0Pct/100) + (Qdis/Q/2)*(cos(2*pi*f0*time)-1); % avg soc vector
thetaAvg = theta0 + zavg*(theta100-theta0); % average stiochiometry vector
iapp = Ipk*sin(2*pi*f0*time); % applied current vector

if ~p.DryRun
    % Generate COMSOL model.
    if p.Verbose
        fprintf('Generating COMSOL model for half-cycle simulation... ');
    end
    modelCOMSOL = genFOM(p.cellModel,'DebugFlag',false);
    if p.Verbose
        fprintf('done!\n');
    end

    % Perform simulation in COMSOL.
    simspec.time = time;
    simspec.mag = Ipk;
    simspec.freq = f0;
    simspec.SOC0 = p.soc0Pct;
    simspec.T = p.TdegC;
    simspec.TSHIFT = 0; % no need to shift for initial discontinuity (sine)
    if p.Verbose
        fprintf('Running half-cycle simulation...\n');
    end
    [~,sim] = simFOM(modelCOMSOL,simspec,'InputType','sin',p.OptSimFOM);
    if p.Verbose
        fprintf('done!\n');
    end
end

simData.time = time(:);
simData.iapp = iapp(:);
if ~p.DryRun
    simData.vcell = sim.Vcell(:);
end
simData.zavg = zavg(:);
simData.thetaAvg = thetaAvg(:);
simData.Qdis = Qdis;
simData.Ipk = Ipk;
simData.f0 = f0;
simData.Tend = Tend;
simData.param = p;
simData.origin__ = 'simHalfCycle';

end