function out = getCellParam(cellModel,varargin)
%GETCELLPARAM Fetch value of cell parameters.
%
% -- Changelog --
% 2023.06.23 | Created | Wesley Hileman <whileman@uccs.edu>

isCellModel = @(x)isstruct(x)&&strcmp(x.type__,'cellModel');
isParamList = @(x)isempty(x)||ischar(x)||iscell(x)||isstruct(x);
parser = inputParser;
parser.addRequired('cellModel',isCellModel);
parser.addOptional('paramList',[],isParamList);
parser.addParameter('TdegC',25,@(x)isscalar(x));
parser.parse(cellModel,varargin{:});
arg = parser.Results;  % structure of validated arguments

% Normalize paramList argument.
paramsStruct = struct;
if isempty(arg.paramList)
    % Get all parameters.
    paramsStruct = cellModel.parameters;
elseif ischar(arg.paramList) || iscell(arg.paramList)
    % Get one for more parameters in list.
    if ~iscell(arg.paramList)
        arg.paramList = strsplit(arg.paramList,{'\\s',','});
    end
    for k = 1:length(arg.paramList)
        paramString = arg.paramList{k};
        parts = strsplit(paramString,'.');
        secName = strtrim(parts{1});
        paramName = strtrim(parts{2});
        if ~isfield(cellModel.parameters,secName)
            error('Section not found in cell model: %s',secName);
        end
        if strcmp(paramName,'*')
            % All parameters in the section.
             paramsStruct.(secName) = cellModel.parameters.(secName);
        elseif ~isfield(cellModel.parameters.(secName),paramName)
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

% Fetch constants needed for Arrhenius temperature correction.
R = TB.const.R;
Tref = cellModel.parameters.const.Tref;
T = arg.TdegC+273.15;

secNames = fieldnames(paramsStruct);
secCount = length(secNames);
for s = 1:secCount
    secName = secNames{s};
    secStruct = paramsStruct.(secName);
    paramNames = fieldnames(secStruct);
    paramCount = length(paramNames);
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
        out.(secName).(paramName) = valueAtT;
    end % for param
end % for sec

if secCount == 1
    if paramCount == 1
        % Single output value.
        out = out.(secNames{1}).(paramNames{1});
    else
        % Single output section.
        out = out.(secNames{1});
    end
end

end