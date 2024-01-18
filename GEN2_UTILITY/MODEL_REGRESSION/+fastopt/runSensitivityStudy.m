function data = runSensitivityStudy(spec,fcn,varargin)
%RUNSENSITIVITYSTUDY Evalulate function sensitivity to parameter values.
%
% data = runSensitivityStudy(spec,fcn) evalulates the MATLAB function
%   FCN for several sets of parameter values specified by SPEC, a structure
%   with the following fields:
%     spec.defaults = structure of unperterbed parameter values
%     spec.singl.values = structure of values to sweep through for each parameter
%       in column-vector format (MSMR parameters can have more than one column)
%     spec.joint.values = (optional) structure specifying joint groups of parameter 
%       values to sweep through. e.g.:
%         joint.grp1.psi = [values for psi]
%         joint.grp1.kD =  [values for kD]
%       Each parameter specified in the joint group must have the same
%       number of values.
%   FCN is a MATLAB function that accepts a structure of parameters as the
%   input and returns a structure of outputs.
%   The return value, DATA, is a structure array with the following fields:
%     data(k).paramname = name of the parameter being perturbed
%     data(k).paramnameEscaped = param name suitable for use in file names
%     data(k).values = values of the parameter evalulated in sens. study
%     data(k).valueLabels = string labels for each parameter value
%     data(k).output = structure array of FCN output for each param value

parser = inputParser;
parser.addRequired('spec',@isstruct);
parser.addRequired('fcn',@(x)isa(x,'function_handle'));
parser.addParameter('Verbose',true,@(x)isscalar(x)&&islogical(x));
parser.parse(spec,fcn,varargin{:});
p = parser.Results; % structure of validated arguments

defaultStruct = fastopt.flattenstruct(spec.defaults);
if isfield(spec, 'singl')
    singlStruct = spec.singl;
else
    singlStruct = struct;
end
if isfield(spec,'joint')
    jointStruct = spec.joint;
else
    jointStruct = struct;
end
singlPerturbTypes = fieldnames(singlStruct);
singlCount = sum(cellfun( ...
    @(x)length(fieldnames(fastopt.flattenstruct(singlStruct.(x)))), ...
        singlPerturbTypes));
jointPerturbTypes = fieldnames(jointStruct);
jointCount = sum(cellfun( ...
    @(x)length(fieldnames(jointStruct.(x))), ...
        jointPerturbTypes));
data.singlCount = singlCount;
data.jointCount = jointCount;
data.totalCount = singlCount+jointCount;

% Evalulate default/baseline output of fcn.
if p.Verbose
    fprintf('Running baseline... ');
end
data.baseline = fcn(spec.defaults);
if p.Verbose
    fprintf('done!\n');
end

getValueFrom.multiplier = @(name, value)defaultStruct.(name)(:)'.*value;
getValueFrom.values = @(name, value)value;
cursor = singlCount+jointCount;

for indPT = 1:length(singlPerturbTypes)
    perturbType = singlPerturbTypes{indPT};
    valueStruct = fastopt.flattenstruct(singlStruct.(perturbType));
    paramNames = fieldnames(valueStruct);
    paramCount = length(paramNames);

    for indJoint = 1:paramCount
        paramName = paramNames{indJoint};
        values = valueStruct.(paramName);

        if p.Verbose
            fprintf('Perturbing %s ', strrep(paramName,'__','.'));
        end
    
        % Evalulate function output for each parameter value.
        clear output; % prevent from persisting between iterations!
        for indVal = length(values):-1:1
            param = defaultStruct;
            param.(paramName) = getValueFrom.(perturbType)(...
                paramName, values(indVal,:));
            param = fastopt.unflattenstruct(param);  % fcn expects hierarchy!
            output(indVal) = fcn(param);
            if p.Verbose
                fprintf('.');
            end
        end % for
    
        if p.Verbose
            fprintf(' done!\n');
        end

        paramNameParts = strsplit(paramName,'__');
        paramNameBase = paramNameParts{end};
    
        data.results(cursor).paramname = strrep(paramName,'__','.');
        data.results(cursor).basename = paramNameBase;
        data.results(cursor).paramnameEscaped = paramName;
        data.results(cursor).values = values;
        if strcmp(perturbType,'multiplier')
            lab = arrayfun(@(x)sprintf('$\\times%.2f$',x),values, ...
                'UniformOutput',false);
        else
            lab = arrayfun(@(x)sprintf('%.2f',x),values, ...
                'UniformOutput',false);
        end
        [n,m] = size(lab);
        if n>1 && m>1
            lab = join(lab,', ');
        end
        data.results(cursor).valuelabels = lab(:);
        data.results(cursor).output = output(:);
        data.results(cursor).perturbType = perturbType;
        data.results(cursor).analysisType = 'singl';
        cursor = cursor - 1;
    end
end

for indPT = 1:length(jointPerturbTypes)
    perturbType = jointPerturbTypes{indPT};
    valuesStruct = jointStruct.(perturbType);
    jointNames = fieldnames(valuesStruct);
    jointCount = length(jointNames);

    for indJoint = 1:jointCount
        jointName = jointNames{indJoint};
        valueStruct = fastopt.flattenstruct(valuesStruct.(jointName));
        paramNames = fieldnames(valueStruct);
        valueCount = length(valueStruct.(paramNames{1}));

        if p.Verbose
            fprintf('Perturbing joint:%s ', jointName);
        end
    
        % Evalulate function output for each set of parameter values.
        clear output; % prevent from persisting between iterations!
        for indVal = 1:valueCount
            param = defaultStruct;
            for indParam = 1:length(paramNames)
                paramName = paramNames{indParam};
                param.(paramName) = getValueFrom.(perturbType)(...
                    paramName, valueStruct.(paramName)(indVal,:));
            end
            param = fastopt.unflattenstruct(param);  % fcn expects hierarchy!
            output(indVal) = fcn(param);
            if p.Verbose
                fprintf('.');
            end
        end % for
    
        if p.Verbose
            fprintf(' done!\n');
        end
    
        valueLabels = cell(valueCount,1);
        for indVal = 1:valueCount
            lab = '';
            for indParam = 1:length(paramNames)
                paramName = paramNames{indParam};
                value = valueStruct.(paramName)(indVal,:);
                if strcmp(perturbType,'multiplier')
                    valStr = sprintf('$\\times%.2f$ ', value);
                else
                    valStr = sprintf('%.2f ', value);
                end
                lab = sprintf('%s, %s', lab, valStr);
            end
            lab = lab(2:end); % remove leading comma
            valueLabels{indVal} = lab;
        end

        basename = '';
        for indParam = 1:length(paramNames)
            paramName = paramNames{indParam};
            parts = strsplit(paramName,'__');
            basename = sprintf('%s, %s', basename, parts{end});
        end
        basename = basename(2:end); % remove leading comma

        data.results(cursor).paramname = ...
            strjoin(strrep(paramNames,'__','.'),',');
        data.results(cursor).basename = basename;
        data.results(cursor).paramnameEscaped = ['joint-' jointName];
        data.results(cursor).values = valueStruct;
        data.results(cursor).valuelabels = valueLabels;
        data.results(cursor).output = output(:);
        data.results(cursor).perturbType = perturbType;
        data.results(cursor).analysisType = 'joint';
        cursor = cursor - 1;
    end
end

end