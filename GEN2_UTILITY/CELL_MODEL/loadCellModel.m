function modelStruct = loadCellModel(filename)
%LOADCELLMODEL Read cell model from Excel spreadsheet.
%
% modelStruct = loadCellModel(filename) reads model parameter values from
%   the spreadsheet MODELNAME and stores them in the structure MODELSTRUCT.
%
% -- Changelog --
% 2023.07.24 | Metadata for diffusion/kinetics model fmt | Wesley Hileman
% 2023.06.21 | Created | Wesley Hileman <whileman@uccs.edu>

parser = inputParser;
parser.addRequired('filename');
parser.parse(filename);
arg = parser.Results;  % structure of validated arguments

[path,name,ext] = fileparts(arg.filename);
if ~strcmpi(ext,'XLSX')
    ext = 'xlsx';
end
filename = fullfile(path,[name '.' ext]);
if ~exist(filename,'file')
    error('Failed to locate cell spreadsheet %s',filename);
end

paramTable = loadSpreadsheet(filename,'Parameters');
sectionStruct = sectionParameterTable(paramTable,arg);
modelStruct = buildModelStruct(sectionStruct,arg);
modelStruct = convertParameterValues(modelStruct,arg);
modelStruct.type__ = 'cellModel';
modelStruct.origin__ = 'loadCellModel';
modelStruct.arg__ = arg;
modelStruct = validateAndPostProcess(modelStruct);

end % loadCellModel()

function dataTable = loadSpreadsheet(filepath,sheetName)
%LOADSPREADSHEET Load contents of Excel spreadsheet into table.

% Read data from Excel sheet into table.
try
    % 'Format','auto' converts everything to text (much faster), only
    % needed for MATLAB >= R2020a.
    dataTable = readtable(filepath, ...
        'Sheet',sheetName,'ReadVariableNames',true,'Format','auto');
catch ME
    if ~strcmp(ME.identifier,'MATLAB:table:parseArgs:BadParamName')
        rethrow(ME);
    end
    % 'Format' not an option in older versions of MATLAB; default
    % behavior is to read everything as text.
    dataTable = readtable(filepath, ...
        'Sheet',sheetName,'ReadVariableNames',true);
end % catch
end % loadSpreadsheet()

function sectionStruct = sectionParameterTable(paramTable,~)
%SECTIONPARAMETERTABLE Divide parameter table into sections.

secioningCol = paramTable.Sectioning;
sectionStruct = struct;

secName = '';
secType = '';
indStart = 0;
for cursor = 1:length(secioningCol)
    item = secioningCol{cursor};
    if ~startsWith(item,'#')
        continue;
    end
    if strcmpi(item,'#end')
        sectionStruct.(secName).name = secName;
        sectionStruct.(secName).type = secType;
        sectionStruct.(secName).paramTable = ...
            paramTable(indStart+1:cursor-1,:);
    else
        parts = strsplit(item(2:end-1),'[');
        secName = parts{1};
        secType = parts{2};
        indStart = cursor;
    end
end

end % sectionParameterTable()

function modelStruct = buildModelStruct(sectionStruct,~)
%BUILDMODEL Create cell model structure from section structure.

% Collect cell/section metadata and parameter properties.
modelStruct = struct;
sectionNames = fieldnames(sectionStruct);
for s = 1:length(sectionNames)
    sec = sectionStruct.(sectionNames{s});
    if strcmpi(sec.type,'Meta')
        % Cell metadata section.
        for k = 1:size(sec.paramTable,1)
            paramName = sec.paramTable.Name{k};
            value = sec.paramTable.Value{k};
            modelStruct.metadata.cell.(paramName) = value;
        end % for
    else
        % Cell parameter section.
        % Names of the parameter properties to retain (exlude the 
        % 'Sectioning' property, which is used only for parsing the file).
        paramPropertyNames = setdiff( ...
            sec.paramTable.Properties.VariableNames,'Sectioning');
        for k = 1:size(sec.paramTable,1)
            paramName = sec.paramTable.Name{k};
            modelStruct.parameters.(sec.name).(paramName) = ... 
                table2struct(sec.paramTable(k,paramPropertyNames));
        end % for

        % Assign metadata for this section.
        modelStruct.metadata.section.(sec.name).name = sec.name;
        modelStruct.metadata.section.(sec.name).type = sec.type;
    end % else
end % for

end % buildModelStruct()

function modelStruct = convertParameterValues(modelStruct, arg)
%CONVERTPARAMETERVALUES Convert parameter-value strings to native MATLAB.

FUNCTION_PARAMETER_NAMES = {'thetas','thetae'};

sectionNames = fieldnames(modelStruct.parameters);
for s = 1:length(sectionNames)
    secName = sectionNames{s};
    sec = modelStruct.parameters.(secName);
    paramNames = fieldnames(sec);
    for p = 1:length(paramNames)
        paramName = paramNames{p};
        param = sec.(paramName);

        % Convert value.
        valueStruct = struct;
        valueString = param.Value;
        valueDouble = str2double(valueString);
        if startsWith(valueString,'#')
            % Lookup table (LUT).
            parts = split(valueString(2:end),'(');
            sheet = parts{1};
            varname = parts{2}(1:end-1); % name of independent var
            tabData = readtable(arg.filename, ...
                'Sheet',sheet, ...
                'ReadVariableNames',false);
            tabData = table2array(tabData);
            tabData1 = tabData(:,1);  % independent variable
            tabData2 = tabData(:,2);  % dependent variable
            valueStruct.type = 'LUT';
            valueStruct.isScalar = false;
            valueStruct.independentVariable = varname;
            valueStruct.x = tabData1;
            valueStruct.y = tabData2;
        elseif ~isnan(valueDouble)
            % Scalar numeric.
            value = valueDouble;
            valueStruct.type = 'Numeric';
            valueStruct.isScalar = true;
            valueStruct.value = value;
        elseif startsWith(valueString,'[') && endsWith(valueString,']')
            % Vector numeric.
            value = eval(valueString);
            value = value(:);
            valueStruct.type = 'Numeric';
            valueStruct.isScalar = false;
            valueStruct.value = value;
        elseif contains(valueString,FUNCTION_PARAMETER_NAMES)
            % Function.
            error(['Not implemented: %s.%s\n' ...
                'Function value not yet implemented.'], ...
                secName,paramName);
        else
            if isempty(valueString)
                valueString = '(empty)';
            end
            error(['Auto-type detection failed: %s.%s\n' ...
                'Value string: %s\n' ...
                'Failed to recognize value as numeric, LUT, or function.'], ...
                secName,paramName,valueString);
        end
        modelStruct.parameters.(secName).(paramName).Value = valueStruct;

        % Convert activation energy.
        EactString = param.Eact;
        EactDouble = str2double(EactString);
        if ~isnan(EactDouble)
            % Scalar numeric.
            Eact = EactDouble;
        elseif startsWith(EactString,'[') && endsWith(EactString,']')
            % Vector numeric.
            Eact = eval(EactString);
            Eact = Eact(:);
        elseif isempty(EactString)
            Eact = 0;
        else
            error(['Auto-type detection failed: %s.%s\n' ...
                'Eact string: %s\n' ...
                'Failed to recognize Eact as numeric scalar or vector.'], ...
                secName,paramName,EactString);
        end
        modelStruct.parameters.(secName).(paramName).Eact = Eact*1000; % convert kJ to J
    end % for
end % for

end % convertParameterValues

function modelStruct = validateAndPostProcess(modelStruct)
%POSTPROCESSANDVALIDATE Validate and post-process cell model.

% Assign chemistry.
modelStruct.metadata.cell.chemistry = 'LMB';

% Determine if this is a standard or lumped model.
if isfield(modelStruct.parameters.const,'W')
    % Lumped-parameter Warburg-resistance model (WRM).
    modelStruct.metadata.cell.type = 'WRM';
    modelStruct.metadata.cell.lumpedParams = true;
else
    % Standard-parameter pseudo-two-dimensional model.
    modelStruct.metadata.cell.type = 'P2DM';
    modelStruct.metadata.cell.lumpedParams = false;
end

% Determine if dll/sep have been combined into single eff layer.
if isfield(modelStruct.parameters,'eff')
    modelStruct.metadata.cell.type = ['RL' modelStruct.metadata.cell.type];
end

% Section post-processing.
sectionNames = fieldnames(modelStruct.metadata.section);
for s = 1:length(sectionNames)
    sec = modelStruct.metadata.section.(sectionNames{s});

    % Additional processing for porous electrode.
    if strcmpi(sec.type,'Electrode3D')
        % Validate OCP / collect OCP metadata.
        isMSMR = false;
        J = 0;
        if all(isfield(modelStruct.parameters.(sec.name),{'U0','X','omega'}))
            isMSMR = true;
            lenU = length(modelStruct.parameters.(sec.name).U0);
            lenX = length(modelStruct.parameters.(sec.name).X);
            lenW = length(modelStruct.parameters.(sec.name).omega);
            if ~(lenU == lenX && lenX == lenW)
                error( ...
                    ['OCP improperly configured: ''%s'' (%s).\n' ...
                     'MSMR parameters U0,X,omega must be of the same ' ...
                     'length.'], ...
                    sec.name,sec.type);
            end
            J = lenU;
            modelStruct.metadata.section.(sec.name).ocp.type = 'msmr';
            modelStruct.metadata.section.(sec.name).ocp.J = J;

            % Ensure lithiation bounds are consistent with Vmin/Vmax definitions!
            if all(isfield(modelStruct.parameters.const,{'Vmin','Vmax'}))
                [vmin, vmax] = getCellParams( ...
                    modelStruct,'const.Vmin,const.Vmax','Output','list');
                [theta100, theta0] = getCellParams( ...
                    modelStruct,'pos.theta100,pos.theta0','Output','list');
                electrodeParams = getCellParams( ...
                    modelStruct,[sec.name '.*']);
                ocpData = MSMR(electrodeParams).ocp('voltage',[vmin vmax]);
                thetamin = ocpData.theta(2);
                thetamax = ocpData.theta(1);
                if theta100<=theta0
                    theta100 = thetamin;
                    theta0 = thetamax;
                else
                    theta100 = thetamax;
                    theta0 = thetamin;
                end
                newParams = struct;
                newParams.(sec.name).theta100 = theta100;
                newParams.(sec.name).theta0 = theta0;
                modelStruct = setCellParam(modelStruct,newParams);
            end
        elseif all(isfield(modelStruct.parameters.(sec.name),{'Uocp','dUocp'}))
            modelStruct.metadata.section.(sec.name).ocp.type = 'explicit';
        else
            error( ...
                ['OCP improperly configured: ''%s'' (%s).\n' ...
                 'Either specify MSMR parameters U0,X,omega OR ' ...
                 'specify Uocp,dUocp explicitly as functions or LUT.'], ...
                sec.name,sec.type);
        end % if

        % Validate kinetics / collect kinetics metadata.
        if isfield(modelStruct.parameters.(sec.name),'k0')
            if isMSMR && length(modelStruct.parameters.(sec.name).k0) == J
                modelStruct.metadata.section.(sec.name).kinetics.type = 'msmr';
            elseif length(modelStruct.parameters.(sec.name).k0) == 1
                modelStruct.metadata.section.(sec.name).kinetics.type = 'bv';
            else
                error( ...
                    ['Kinetics improperly configured: ''%s'' (%s).\n' ...
                     'For MSMR models, k0 should be a vector of length J. ' ...
                     'Otherwise k0 should be a scalar.'], ...
                    sec.name,sec.type);
            end
        elseif all(isfield(modelStruct.parameters.(sec.name),{'k0Spline','k0Theta'}))
            modelStruct.metadata.section.(sec.name).kinetics = 'spline';
        elseif all(isfield(modelStruct.parameters.(sec.name),{'k0Linear','k0Theta'}))
            modelStruct.metadata.section.(sec.name).kinetics = 'linear';
        else
            error( ...
                ['Kinetics improperly configured: ''%s'' (%s).\n' ...
                 'Specify the lumped reaction rate-constant k0 OR ' ...
                 'k0Spline and k0Theta OR k0Linear and k0Theta.'], ...
                sec.name,sec.type);
        end % if

        % Validate diffusivity / collect diffusivity metadata.
        if isfield(modelStruct.parameters.(sec.name),'Dsref')
            % OK: can still use Baker-Vergrugge (MSMR) diffusivity model
            % if ~isMSMR
            %     error( ...
            %         ['Diffusivity improperly configured: ''%s'' (%s).\n' ...
            %          'Non-MSMR models cannot use Dsref; specify Ds instead.'], ...
            %         sec.name,sec.type);
            % end
             modelStruct.metadata.section.(sec.name).solidDiffusion.type = 'msmr';
        elseif isfield(modelStruct.parameters.(sec.name),'Ds')
            modelStruct.metadata.section.(sec.name).solidDiffusion.type = 'fixed';
        elseif all(isfield(modelStruct.parameters.(sec.name),{'DsSpline','DsTheta'}))
            modelStruct.metadata.section.(sec.name).solidDiffusion.type = 'spline';
        elseif all(isfield(modelStruct.parameters.(sec.name),{'DsLinear','DsTheta'}))
            modelStruct.metadata.section.(sec.name).solidDiffusion.type = 'linear';
        else
            error( ...
                ['Diffusivity improperly configured: ''%s'' (%s).\n' ...
                 'Specify the diffusivity coefficient Dsref OR ' ...
                 'fixed diffusivity Ds OR DsSpline and DsTheta OR ' ...
                 'DsLinear and DsTheta.'], ...
                sec.name,sec.type);
        end % if
    end % if Electrode3D
end % for
end % validateAndPostProcess()