function simData = simEIS(varargin)
%SIMEIS Run full-order EIS simulation in COMSOL over specified frequencies.
%
% -- Usage --
% data = SIMEIS(cellModel,freq,socPct) simulates the full-order electrochemical
%   impedance spectroscopy (EIS) response a cell in COMSOL over the 
%   frequencies provided in the vector FREQ. CELLMODEL is a model of the 
%   cell obtained from LOADCELLPARAMS. SOCPCT is the state-of-charge of the
%   cell in percent [%]; if SOCPCT is a vector, impedance spectra is
%   simulayed at multiple SOC setpoints. DATA is a structure containing 
%   the time-domain simulation data (see description below).
%
% data = SIMEIS(cellModel,freq,socPct,TdegC) also specifies the temperature
%   at which to run the simulation. Default TdegC=25.
%
% data = SIMEIS(...,'Vars',varsStruct) specifies the electrochemical variables 
%   to save in addition to Vcell and Iapp as well as the x-locations at which
%   to evalulate those variables. varsStruct is a structure containing
%   fields for each variable. The value of a field specifies the
%   x-locations at which to evalulate the specified variable.
%   Use the string value of 'mesh' to save all available x-locations.
%   Default varsStruct=struct('Phise',[0 3],'Phie',3).
%   NOTE: this function does not interpolate; if an x-location is specified
%   that is not in the COMSOL mesh, the closest x-location to that provided
%   is fetched (see which x-locations were fetched by accessing the `xlocs`
%   field of the returned DATA structure).
%
% data = SIMEIS(...,'I',I) specifies the amplitude of the input sinusoid 
%   as I [A]. Default I=C/10.
%
% data = SIMEIS(...,'Ns',Ns) specifies the sampling rate as the number of 
%   samples per period of the sinusoid [Sa/period]. Default Ns=40.
%
% data = SIMEIS(...,'Nt',Nt) specifies the length of the transient interval 
%   in periods of the applied sinusoid [periods]. Default Nt=20.
%
% data = SIMEIS(...,'Nss',Nss) specifies the length of the steady-state 
%   interval in periods of the applied sinusoid [periods]. Default Nss=10.
%
% data = SIMEIS(...,'Verbose',true) outputs status messages to the command
%   window.
%
% The output, DATA, is a structure with the following fields:
%    ss = structure array containing steady-state waveforms. ss(kf,kz) is
%      a structure with the following fields corresponding to the kf-th
%      frequency in the input frequency vector freq and the kz-th SOC in
%      the input SOC setpoint vector socPct:
%      - time: simulation time
%      - Vcell: vector of cell voltage vs time
%      - Iapp: vector of applied current vs time
%      - {Vars}: fields for each variable specified in the 'Vars' input
%        cell array. The value is a matrix whose first dimension (rows)
%        correspond to time and whose second dimension corresponds to
%        x-locations (columns).
%    xlocs = structure mapping variables to x-locations in the cell
%      sandwich.
%    param = structure of parameters supplied to the function.
%
% -- Changelog --
% 2023.06.11 | Support multiple SOC setpoints | Wesley Hileman
% 2023.04.04 | Created | Wesley Hileman <whileman@uccs.edu>

% Fetch and validate input arguments.
isinteger = @(x)floor(x)==x;
parser = inputParser;
parser.addRequired('cellModel',@(x)isscalar(x)&&isstruct(x));
parser.addRequired('freq',@(x)isvector(x)&&all(x>0));
parser.addRequired('socPct',@(x)isvector(x)&&all(0<=x&x<=100));
parser.addOptional('TdegC',25,@(x)isscalar(x));
parser.addParameter('Vars',struct('Phise',[0 3],'Phie',3),@(x)isstruct(x));
parser.addParameter('I',[],@(x)isscalar(x)&&x>0);
parser.addParameter('Ns',40,@(x)isscalar(x)&&isinteger(x)&&x>0);
parser.addParameter('Nt',20,@(x)isscalar(x)&&isinteger(x)&&x>0);
parser.addParameter('Nss',10,@(x)isscalar(x)&&isinteger(x)&&x>0);
parser.addParameter('Verbose',false,@(x)isscalar(x)&&islogical(x));
parser.addParameter('OptSimFOM',struct,@isstruct);
parser.parse(varargin{:});
p = parser.Results;  % structure of validated arguments

if isfield(p.OptSimFOM,'VcellOnly') && p.OptSimFOM.VcellOnly
    p.Vars = struct; % only cell voltage is available!
end

if isempty(p.I)
    % Use C/10 rate as default.
    p.I = getCellParams(p.cellModel,'const.Q')/10;
end

% Generate COMSOL model.
if p.Verbose
    fprintf('Generating COMSOL model for EIS simulation...\n');
end
modelCOMSOL = genFOM(p.cellModel);

% Run simulation at each SOC setpoint and frequency (backwards to avoid 
% need to pre-allocate structure arrays).
ff = p.freq(:);
varnames = fieldnames(p.Vars);
clear ssdata;
xlocsStruct = struct;
for kz = length(p.socPct):-1:1
    for kf = length(ff):-1:1
        f0 = ff(kf);
    
        % Generate time vector.
        kk = [0:p.Nt*p.Ns,(p.Nt*p.Ns+1):(p.Nt+p.Nss)*p.Ns]; % Discrete-time vector [sample number].      
        fs = p.Ns*f0; % Sampling rate [Sa/s].
        tt = kk/fs; % Time vector [s].
        ss = kk>p.Nt*p.Ns; % Logical indicies to steady-state interval.
    
        % Perform simulation in COMSOL.
        simspec.time = tt;
        simspec.mag = p.I;
        simspec.freq = f0;
        simspec.SOC0 = p.socPct(kz);
        simspec.T = p.TdegC;
        simspec.TSHIFT = 0; % no need to shift for initial discontinuity (sine)
        if p.Verbose
            fprintf('Running SOC=%8.3%% f=%8.3gHz (%d of %d)...\n', ...
                p.socPct(kz),f0,length(ff)-kf+1,length(ff));
        end
        [~,sim] = simFOM(modelCOMSOL,simspec,'InputType','sin',p.OptSimFOM);
        if p.Verbose
            fprintf('done!\n');
        end
    
        % Save steady-state waveforms.
        ssdata(kz,kf).param.f0 = f0;
        ssdata(kz,kf).param.fs = fs;
        ssdata(kz,kf).param.N = length(tt(ss));
        ssdata(kz,kf).param.socPct = p.socPct(kz);
        ssdata(kz,kf).time = tt(ss);
        ssdata(kz,kf).time = ssdata(kz,kf).time - ssdata(kz,kf).time(1); % start time at zero
        ssdata(kz,kf).Iapp = sim.Iapp(ss);
        ssdata(kz,kf).Vcell = sim.Vcell(ss);
        for i = 1:length(varnames)
            varname = varnames{i};
            timePositionData = sim.(varname);
            xlocs = sim.xLocs.(varname);      % COMSOL mesh locations
            xlocsDesired = p.Vars.(varname);  % desired x-locations
            if (ischar(xlocsDesired)||isstring(xlocsDesired))&&strcmpi(xlocsDesired,'mesh')
                % Store variable at all mesh locations.
                ssdata(kz,kf).(varname) = timePositionData(ss,:);
            else
                [~,indxx] = min(abs(xlocs(:)-xlocsDesired(:)'));
                ssdata(kz,kf).(varname) = timePositionData(ss,indxx);
                xlocs = xlocs(indxx); % true x-locations
            end
            if ~isfield(xlocsStruct,varname)
                xlocsStruct.(varname) = xlocs;
            end
        end
    end % for freq
end % for soc

if p.Verbose
    fprintf('Finished full-order EIS simulation\n');
end

% Collect results.
simData.ss = ssdata;
simData.xlocs = xlocsStruct;
simData.param = p;
simData.origin__ = 'simEIS';

end