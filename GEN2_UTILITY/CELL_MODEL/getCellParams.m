function varargout = getCellParams(cellModel,varargin)
%GETCELLPARAMS Fetch values of cell parameters at given temperature/SOC.
%
% -- Usage --
% outStruct = GETCELLPARAMS(cellModel) evalulates all parameters of the
%   cell model CELLMODEL at T=25degC and SOC=50% and returns a structure
%   containing the parameter values, OUTSTRUCT.
%
% outStruct = GETCELLPARAMS(cellModel,paramList) evalulates the subset of
%   parameters sepecified by PARAMLIST. PARAMLIST may be a character vector
%   or a cell array of character vectors. See below for character-vector
%   format.
%
% outStruct = GETCELLPARAMS(...,'TdegC',TdegC) evalulates parameter values
%   at temperature TdegC instead of the default 25degC.
%
% outStruct = GETCELLPARAMS(...,'socPct',socPct) evalulates parameter
%   values at SOC socPct instead of the default 50%.
%
% outStruct = GETCELLPARAMS(...,'Output,'nested') forces OUTSTRUCT to
%   appear in the same nested structure format as the cell model 
%   (default is to return a scalar when a single parameter is requested, 
%   flat structure when parameters in a single section are requested, and
%   postfix flat structure when single type of parameter requested from
%   multiple sections).
%
% outStruct = GETCELLPARAMS(...,'Output','postfix') flattens the output
%   structure by post-fixing section names to parameter names.
%
% [out1,out2,out3,...] = GETCELLPARAMS(...,'Output','list') returns the
%   requsted parameter values as a list of output arguments.
%
% Format for paramList character vector:
%    Single parameter: '(sectionName).(parameterName)'
%    Multiple parameters: '(sec1).(param1) (sec2).(param2) (sec3).(param3)'
%    All parameters in section: '(sectionName).*'
%    All sections containing parameter: '*.(parameterName)'
%
% -- Examples --
% Let p2dm be a Gen2 P2DM structure.
% Let llpm be a Gen1 LLPM sructure.
%
% paramStruct = getCellParams(p2dm); % get all params
% lenStruct = getCellParams(p2dm,'*.L'); % get region lengths
% posStruct = getCellParams(p2dm,'pos.*'); % get porous-electrode params
% Q = getCellParams(llpm,'const.Q'); % get single parameter
% paramStruct = getCellParams(llpm,'*.kappa,neg.*,const.Q'); % get multiple
%
% -- Changelog --
% 2023.06.23 | Created | Wesley Hileman <whileman@uccs.edu>

isParamList = @(x)isempty(x)||iscell(x)||isstruct(x)||...
    ischar(x)&&contains(x,'.');
isOutputMode = @(x)any(strcmpi(x,{'auto','nested','postfix','list'}));
parser = inputParser;
parser.addRequired('cellModel',@isstruct);
parser.addOptional('paramList',[],isParamList);
parser.addParameter('TdegC',25,@(x)isscalar(x));
parser.addParameter('socPct',50,@(x)isscalar(x)&&0<=x&&x<=100);
parser.addParameter('Output','auto',isOutputMode);
parser.parse(cellModel,varargin{:});
arg = parser.Results;  % structure of validated arguments

% Determine type of model input.
cellModelType = getCellModelType(arg.cellModel,false);

% Normalize paramList argument.
paramsStruct = struct;
if isempty(cellModelType)
    % Structure of parameter values.
    cellParams = arg.cellModel;
elseif strcmp(cellModelType,'LLPM')
    % Gen1 LPM.
    cellParams = arg.cellModel.function;
else
    % Gen2 model.
    cellParams = arg.cellModel.parameters;
end
if isempty(arg.paramList)
    % Get all parameters.
    paramsStruct = cellParams;
elseif ischar(arg.paramList) || iscell(arg.paramList)
    % Get one or more parameters in list.
    if ~iscell(arg.paramList)
        arg.paramList = strsplit( ...
            arg.paramList,{'\s',',',';'}, ...
            'DelimiterType', 'RegularExpression');
    end
    for k = 1:length(arg.paramList)
        paramString = arg.paramList{k};
        parts = strsplit(paramString,'.');
        secName = strtrim(parts{1});
        paramName = strtrim(parts{2});
        if strcmp(secName,'*')
            % All sections containing the parameter.
            secNames = fieldnames(cellParams);
            for s = 1:length(secNames)
                secName = secNames{s};
                secStruct = cellParams.(secName);
                if isfield(secStruct,paramName)
                    paramsStruct.(secName).(paramName) = true;
                end % if
            end % for sec
        elseif ~isfield(cellParams,secName)
            error('Section not found in cell model: %s',secName);
        elseif strcmp(paramName,'*')
            % All parameters in the section.
            paramsStruct.(secName) = cellParams.(secName);
        elseif ~isfield(cellParams.(secName),paramName)
            error('Parameter not found in cell model: %s.%s', ...
                secName,paramName);
        else
            % Single parameter in the section.
            paramsStruct.(secName).(paramName) = true;
        end % if
    end % for paramList
elseif isstruct(arg.paramList)
    paramsStruct = arg.paramList;
end

if isempty(cellModelType)
    % Structure of parameter values.
    [outStruct, totalParamCount] = getStructParams(arg,paramsStruct);
elseif strcmp(cellModelType,'LLPM')
    % Gen1 legacy LPM.
    [outStruct, totalParamCount] = getGen1Params(arg,paramsStruct);
else
    % Gen2 model.
    [outStruct, totalParamCount] = getGen2Params(arg,paramsStruct);
end

% Format output.
if strcmpi(arg.Output,'auto')
    secNames = fieldnames(outStruct);
    secCount = length(secNames);
    if secCount == 1
        paramNames = fieldnames(outStruct.(secNames{1}));
        paramCount = length(paramNames);
        if paramCount == 1
            % Single output value.
            varargout{1} = outStruct.(secNames{1}).(paramNames{1});
        else
            % Single output section.
            varargout{1} = outStruct.(secNames{1});
        end % if
    else
        % Multiple output sections.
        sectionsHaveSingleParam = true;
        sectionParamNames = cell(1,secCount);
        for s = 1:secCount
            secName = secNames{s};
            secStruct = outStruct.(secName);
            paramNames = fieldnames(secStruct);
            paramCount = length(paramNames);
            if paramCount ~= 1
                sectionsHaveSingleParam = false;
                break;
            end
            sectionParamNames{s} = paramNames{1};
        end % for section
        if sectionsHaveSingleParam && isequal(sectionParamNames{:})
            % Same parameter from each section; flat postfix output.
            varargout{1} = postfixFlatten(outStruct);
        else
            % Mixed parameters from each; structure output.
            varargout{1} = outStruct;
        end % if
    end % if
elseif strcmpi(arg.Output,'nested')
    % Structure output.
    varargout{1} = outStruct;
elseif strcmpi(arg.Output,'postfix')
    % Flat structure output (fieldnames postfixed with section names).
    varargout{1} = postfixFlatten(outStruct);
elseif strcmpi(arg.Output,'list')
    % Multiple output arguments as list.
    cursor = 1;
    varargout = cell(1,totalParamCount);
    secNames = fieldnames(outStruct);
    secCount = length(secNames);
    for s = 1:secCount
        secName = secNames{s};
        secStruct = paramsStruct.(secName);
        paramNames = fieldnames(secStruct);
        paramCount = length(paramNames);
        for p = 1:paramCount
            paramName = paramNames{p};
            varargout{cursor} = outStruct.(secName).(paramName);
            cursor = cursor + 1;
        end % for param
    end % fo sec
end % if

end

function [outStruct, totalParamCount] = getStructParams(arg,paramsStruct)

totalParamCount = 0;
secNames = fieldnames(paramsStruct);
secCount = length(secNames);
for s = 1:secCount
    secName = secNames{s};
    secStruct = paramsStruct.(secName);
    paramNames = fieldnames(secStruct);
    paramCount = length(paramNames);
    totalParamCount = totalParamCount + paramCount;
    for p = 1:paramCount
        paramName = paramNames{p};
        paramValue = arg.cellModel.(secName).(paramName);
        outStruct.(secName).(paramName) = paramValue;
    end % for param
end % for sec

end

function [outStruct, totalParamCount] = getGen1Params(arg,paramsStruct)

T = arg.TdegC+273.15;
thetap = arg.cellModel.function.pos.soc(arg.socPct/100);
totalParamCount = 0;

secNames = fieldnames(paramsStruct);
secCount = length(secNames);
for s = 1:secCount
    secName = secNames{s};
    secStruct = paramsStruct.(secName);
    paramNames = fieldnames(secStruct);
    paramCount = length(paramNames);
    totalParamCount = totalParamCount + paramCount;
    for p = 1:paramCount
        paramName = paramNames{p};
        paramFcn = arg.cellModel.function.(secName).(paramName);
        outStruct.(secName).(paramName) = paramFcn(thetap,T);
    end % for param
end % for sec

end

function [outStruct, totalParamCount] = getGen2Params(arg,paramsStruct)

% Fetch constants needed for Arrhenius temperature correction.
R = TB.const.R;
Tref = arg.cellModel.parameters.const.Tref.Value.value;
T = arg.TdegC+273.15;
totalParamCount = 0;

% Fetch parameters.
secNames = fieldnames(paramsStruct);
secCount = length(secNames);
for s = 1:secCount
    secName = secNames{s};
    secStruct = paramsStruct.(secName);
    paramNames = fieldnames(secStruct);
    paramCount = length(paramNames);
    totalParamCount = totalParamCount + paramCount;
    for p = 1:paramCount
        paramName = paramNames{p};
        param = arg.cellModel.parameters.(secName).(paramName);

        % Decode value.
        if strcmp(param.Value.type,'Numeric')
            valueAtTref = param.Value.value;
        else
            error( ...
                ['Not implemented: %s.%s\n' ...
                 'Decoder for value type ''%s'' not implemented.'], ...
                 secName,paramName,paramsStruct.Value.type);
        end % if

        % Apply Arrhenius temperature correction (if needed).
        Eact = param.Eact;
        if Eact == 0
            valueAtT = valueAtTref;
        else
            valueAtT = valueAtTref.*exp(Eact*(1/Tref-1/T)/R);
        end % if

        % Store value in output structure.
        outStruct.(secName).(paramName) = valueAtT;
    end % for param
end % for sec

end

function postfixStruct = postfixFlatten(nestedStruct)

postfixStruct = struct;
secNames = fieldnames(nestedStruct);
secCount = length(secNames);
for s = 1:secCount
    secName = secNames{s};
    secStruct = nestedStruct.(secName);
    paramNames = fieldnames(secStruct);
    paramCount = length(paramNames);
    for p = 1:paramCount
        paramName = paramNames{p};
        if strcmp(secName,'const')
            fieldName = paramName;
        else
            fieldName = [paramName secName];
        end
        postfixStruct.(fieldName) = nestedStruct.(secName).(paramName);
    end % for param
end % for sec

end