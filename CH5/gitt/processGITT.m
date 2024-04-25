function out = processGITT(gittIn,ocpIn,varargin)
%PROCESSGITT Extract diffusivity estimates from a laboratory GITT sequence.
%  Works for half-cells and lithium-metal battery cells with a single
%  porous electrode.
%
%  Based on the relaxation technique presented by:
%    Stephen Dongmin Kang and William C. Chueh 2021 
%    J. Electrochem. Soc. 168 120504
%
% -- Usage --
% out = processGITT(gittIn,ocpIn) segments the specified GITT sequence 
%   into constant-current and relaxation intervals and them performs 
%   linear regression on the relaxation voltage to estimate the solid 
%   diffusivity.
%   - gittIn is a structure of the GITT data. Required fields:
%       gittIn.time     Time vector [s]
%       gittIn.Vcell    Cell voltage vector [V]
%       gittIn.Iapp     Applied current vector, >0 for discharge [A]
%       gittIn.soc0Pct  Starting SOC of the cell [%]
%       gittIn.TdegC    Cell temperature [degC]
%   - ocpIn is a structure of OCP data used to establish the SOC and
%     fractional lithiation of the cell. Required fields:
%       ocpIn.QAh       Capacity of the cell [Ah]
%       ocpIn.theta0    Stiochiometry at 0% SOC [fractional]
%       ocpIn.theta100  Stiochiometry at 100% SOC [fractional]
%       ocpIn.U0        MSMR parameters
%       ocpIn.X         "
%       ocpIn.omega     "
%   - out is a structure of output data with the following fields:
%       out.segments    Structure array of processed GITT segments
%       out.arg         Structure of arguments supplied to this function.
%
% out = processGITT(...,'tcPct',tcPct) sets the cutoff time for estimating
%   the solid diffusivity as a percentage of the total relaxation interval.
%   Data for trelax>=(tcPct/100)*Tr, where Tr is the total relaxation time
%   and trelax denotes the time variable since the relaxation began, will
%   be used to estimate Ds. DEFAULT: 30%.
%
% -- Changelog --
% 2023.10.10 | Created | Wesley Hileman <whileman@uccs.edu>

parser = inputParser;
parser.addRequired('gittIn',@(x)isstruct(x)&&isscalar(x));
parser.addRequired('ocpIn',@(x)isstruct(x)&&isscalar(x));
parser.addParameter('tcPct',30,@(x)isnumeric(x)&&isscalar(x)&&0<=x&&x<=100);
parser.parse(gittIn,ocpIn,varargin{:});
arg = parser.Results;

% Collect GITT and OCP data.
time = gittIn.time(:).';
Iapp = gittIn.Iapp(:).';
Vcell = gittIn.Vcell(:).';
TdegC = gittIn.TdegC;
soc0Pct = gittIn.soc0Pct;
QtotAh = ocpIn.QAh;
theta0 = ocpIn.theta0;
theta100 = ocpIn.theta100;

% Segment pulses ----------------------------------------------------------
I = max(abs(Iapp));       % pulse magnitude [A]
relax = abs(Iapp)<I/100;  % logical indicies to relaxation intervals
% For explaination of lines below, let: 
% relax = [ 0  0  0  1  1  1  0  0  0  1  1  1]
% index:    1  2  3  4  5  6  7  8  9 10 11 12
edges = [0 diff(relax)];  % => [ 0  0  0  1  0  0 -1  0  0  1  0  0]
indRelaxStart = find(edges==1);  % => [4 10]
indRelaxEnd = [find(edges==-1)-1 length(edges)];  % => [6 12]
% Create struct array locating the bounds of each segment.
clear segments;
for k = length(indRelaxStart):-1:1
    if k == 1
        indPulseStart = 1; % first pulse starts at time=0
    else
        indPulseStart = indRelaxEnd(k-1)+1;
    end
    indPulseEnd = indRelaxStart(k)-1;
    indPulse = indPulseStart:indPulseEnd;
    indRelax = indRelaxStart(k):indRelaxEnd(k);
    tpulse = time(indPulse);
    tpulse = tpulse - tpulse(1);
    trelax = time(indRelax);
    trelax = trelax - trelax(1);
    tau = tpulse(end) - tpulse(1);
    Iavg = mean(Iapp(indPulse));
    QdisAh = trapz(tpulse,Iapp(indPulse))/3600;
    segments(k).indPulse = indPulse;
    segments(k).indRelax = indRelax;
    segments(k).tau = tau;
    segments(k).trelax = trelax;
    segments(k).Vrelax = Vcell(indRelax);
    segments(k).Iavg = Iavg;
    segments(k).QdisAh = QdisAh;
    segments(k).socPct = [];  % allocate space for later assignment
    segments(k).theta = [];
    segments(k).slope = [];
    segments(k).intercept = [];
    segments(k).DsRel = [];
    segments(k).dUocp = [];
end
% Determine SOC and lithiation at start of each relaxation interval.
socPct = soc0Pct;
for k = 1:length(segments)
    socPct = socPct - 100*segments(k).QdisAh/QtotAh;
    theta = theta0 + (socPct/100)*(theta100 - theta0);
    segments(k).socPct = socPct;
    segments(k).theta = theta;
end

% Process segments --------------------------------------------------------
ocpData = MSMR(ocpIn).ocp('theta',[segments.theta],'TdegC',TdegC);
for k = 1:length(segments)
    seg = segments(k);
    dUocp = ocpData.dUocp(k);
    trelax = seg.trelax;
    tdiff = sqrt(trelax + seg.tau) - sqrt(trelax);
    Vrelax = seg.Vrelax;
    % Isolate linear portion of the response.
    ind = trelax>=(arg.tcPct/100)*trelax(end);
    t = tdiff(ind);
    v = Vrelax(ind);
    % Estimate slope of the linear portion.
    coeffs = [t(:) ones(length(t),1)]\v(:);
    slope = coeffs(1);
    intercept = coeffs(2);
    % Compute relative diffusivity from slope.
    DsRel = (dUocp/slope)^2;
    % Store results.
    segments(k).slope = slope;
    segments(k).intercept = intercept;
    segments(k).DsRel = DsRel;
    segments(k).dUocp = dUocp;
end

out.segments = segments;
out.arg = arg;
end