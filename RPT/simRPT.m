function simData = simRPT(varargin)
%SIMRPT Simulate reference performance test (RPT) in COMSOL.
%
% The RPT consists of a half-cycle (dis)charge and (NL)EIS.
%
% -- Usage --
% simData = simRPT(cellModel, ArgHalfCycle, ArgEIS) simulates an
%   RPT in COMSOL by calling simHalfCycle and simEIS with the
%   respective argument structures, ArgHalfCycle and ArgEIS.
%   The first argument is automatically set to cellModel.
%
% -- Changelog --
% 2023.06.11 | Created | Wesley Hileman <whileman@uccs.edu>

parser = inputParser;
parser.addRequired('cellModel',@isstruct);
parser.addRequired('ArgHalfCycle',@iscell);
parser.addRequired('ArgEIS',@iscell);
parser.addParameter('Verbose',false,@(x)islogical(x)&&isscalar(x));
parser.addParameter('OptSimFOM', ...
    struct('VcellOnly',true,'DebugFlag',false),@isstruct)
parser.parse(varargin{:});
p = parser.Results;  % structure of validated arguments.

% Build structures of keyword arguments to pass to simYYY functions.
OptHalfCycle.OptSimFOM = p.OptSimFOM;
OptEIS.OptSimFOM = p.OptSimFOM;

if p.Verbose
    fprintf('Running RPT\n==============\n');
end
simData.halfCycle = simHalfCycle( ...
    p.cellModel,p.ArgHalfCycle{:},OptHalfCycle);
simData.eis = simEIS( ...
    p.cellModel,p.ArgEIS{:},OptEIS);
simData.param = p;
simData.origin__ = 'simRPT';

end