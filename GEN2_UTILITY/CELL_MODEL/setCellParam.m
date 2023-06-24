function cellModel = setCellParam(cellModel,params)
%SETCELLPARAM Update values of numerical cell parameters.
%
% newModel = setCellParam(oldModel,paramStruct) replaces the values of
%   numerical parameters in the cell model OLDMODEL with the new values
%   specified in the structure PARAMSTRUCT and returns a new cell model
%   NEWMODEL.
%
% If any parameter is not a numeric (e.g., a LUT), an error is thrown. 
%
% -- Changelog --
% 2023.06.23 | Created | Wesley Hileman <whileman@uccs.edu>

parser = inputParser;
parser.addRequired('cellModel',@isstruct);
parser.addRequired('params',@isstruct);
parser.parse(cellModel,params);
arg = parser.Results;  % structure of validated arguments

secNames = fieldnames(arg.params);
for s = 1:length(secNames)
    secName = secNames{s};
    sec = arg.params.(secName);
    paramNames  = fieldnames(sec);
    for p = 1:length(paramNames)
        paramName = paramNames{p};
        newValue = sec.(paramName);
        param = arg.cellModel.parameters.(secName).(paramName);
        if ~strcmp(param.Value.type,'Numeric')
            error('Type error: %s.%s is not numeric!',secName,paramName);
        end
        cellModel.parameters.(secName).(paramName).Value.value = newValue;
    end % for param
end % for section

end