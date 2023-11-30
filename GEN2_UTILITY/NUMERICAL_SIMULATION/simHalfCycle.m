function simData = simHalfCycle(varargin)
%SIMHALFCYCLE Simulate half-cycle dis/charge of a LMB cell in COMSOL. Useful
%  for comparison to the perturbation-resistance approximation developed by 
%  Baker and Verbrugge [1].
%
% -- Usage --
% simData = simHalfCycle(cellModel, IavgC, soc0Pct, socfPct) simulates the
%   response of the full-order cell model CELLMODEL to a sinusoidal
%   half-wave dis/charge current with average value IAVGC [C-rate]. The 
%   initial SOC is set by SOC0PCT [%] and the final SOC by SOCFPCT [%].
%   For negative C-rates (implies charge), SOC0PCT and SOCFPCT are reversed if 
%   SOCFPCT < SOC0PCT. IAVGC may be a vector, in which case a structure
%   array of output data is returned with an element for each Iavg.
%   **The COMSOL LiveLink interface for MATLAB must be running.**
%
% simData = simHalfCycle(cellModel, Iavg, soc0Pct, socfPct, TdegC) also
%   specifies the temperature [degC] at which the simulation is run.
%
% simData = simHalfCycle(...,'Ns',Ns) uses Ns samples instead of the
%   default 1000.
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
parser.addRequired('IavgC',@(x)isvector(x));
parser.addRequired('soc0Pct',@(x)isscalar(x)&&0<=x&&x<=100);
parser.addRequired('socfPct',@(x)isscalar(x)&&0<=x&&x<=100);
parser.addOptional('TdegC',25,@(x)isscalar(x));
parser.addParameter('Ns',1000,@(x)isscalar(x)&&x>0);
parser.addParameter('Verbose',false,@(x)isscalar(x)&&islogical(x));
parser.addParameter('OptSimFOM',struct,@isstruct);
parser.addParameter('DryRun',false,@(x)isscalar(x)&&islogical(x));
parser.parse(varargin{:});
p = parser.Results;  % structure of validated params

cellModelType = getCellModelType(p.cellModelOrQ,false);
if isempty(cellModelType)
    cellModel = [];
    Q = p.cellModelOrQ;  % not a cell model, interpret as capacity!
    theta0 = NaN;
    theta100 = NaN;
else
    cellModel = convertCellModel(p.cellModelOrQ,'LLPM');
    [Q,theta0,theta100] = getCellParams(cellModel, ...
        'const.Q pos.theta0 pos.theta100','Output','list');
end
Iavg = Q*p.IavgC;  % average iapp in Amperes

if ~p.DryRun
    % Generate COMSOL model.
    if p.Verbose
        fprintf('Generating COMSOL model for half-cycle simulation... ');
    end
    modelCOMSOL = genFOM(cellModel,'DebugFlag',false);
    if p.Verbose
        fprintf('done!\n');
    end
end

% Perform simulations at each C-rate.
clear hcycSeries;
for k = length(Iavg):-1:1
    % Ensure SOC limits are correct given sign of iavg.
    iavg = Iavg(k);
    z0 = p.soc0Pct/100;
    zf = p.socfPct/100;
    if iavg<0 && zf<z0
        tmp = z0;
        z0 = zf;
        zf = tmp;
    end
    iavg = abs(iavg); % SOC limits now reflect the sign

    % Compute cell capacity to discharge (will be negative for charge).
    Qdis = Q*(z0-zf);   % capacity to discharge [Ah] 

    % Compute half-cycle waveform parameters.
    Ipk = sign(Qdis)*iavg*pi/2;   % amplitude of iapp half-sinusoid [A]
    f0 = iavg/abs(Qdis)/3600/2;   % cyclic frequency of iapp sinusoid [Hz]
    Tend = 1/f0/2;                % time at end of half-cycle [s]
    time = linspace(0,Tend,p.Ns); % time vector [s]      
    zavg = (p.soc0Pct/100) + (Qdis/Q/2)*(cos(2*pi*f0*time)-1); % avg soc vector
    thetaAvg = theta0 + zavg*(theta100-theta0); % average stiochiometry vector
    iapp = Ipk*sin(2*pi*f0*time); % applied current vector
    
    if ~p.DryRun
        % Perform simulation in COMSOL.
        simspec.time = time;
        simspec.mag = Ipk;
        simspec.freq = f0;
        simspec.SOC0 = z0*100;
        simspec.T = p.TdegC;
        simspec.TSHIFT = 0; % no need to shift for initial discontinuity (sine)
        if p.Verbose
            fprintf('Running half-cycle simulation @ Iavg=%5.2fC (#%d)...',p.IavgC(k),k);
        end
        [~,sim] = simFOM(modelCOMSOL,simspec,'InputType','sin',p.OptSimFOM);
        if p.Verbose
            fprintf('done!\n');
        end
    end % if
    
    hcycSeries(k).time = time(:);
    hcycSeries(k).iapp = iapp(:);
    if ~p.DryRun
        hcycSeries(k).vcell = sim.Vcell(:);
    end
    hcycSeries(k).zavg = zavg(:);
    hcycSeries(k).thetaAvg = thetaAvg(:);
    hcycSeries(k).Qdis = Qdis;
    hcycSeries(k).Ipk = Ipk;
    hcycSeries(k).f0 = f0;
    hcycSeries(k).Tend = Tend;
    hcycSeries(k).z0 = z0;
    hcycSeries(k).zf = zf;
end

simData.series = hcycSeries;
simData.IavgC = p.IavgC;
simData.TdegC = p.TdegC;
simData.cellModel = cellModel;
simData.Q = Q;
simData.arg = p;
simData.origin__ = 'simHalfCycle';

end