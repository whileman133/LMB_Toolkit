function obj = reverseFcn(fcn)
    %REVERSEFCN Determine the parameter value associated with a compiled
    % inline function from loadCellParams().
    %
    % obj = REVERSEFCN(fcn) explodes the inline function FCN and returns an
    %  object OBJ (a structure) containing the useful bits. OBJ has three
    %  fields:
    %    type: 'LUT', 'function', 'vector', or 'scalar'
    %    value: value, format depends on the type
    %    Ea: activation energy (scalar or vector, zero if not given) [J/mol]
    %
    % Three types of parameters are currently supported by readCellParams:
    %   1. Lookup tables (LUTs)
    %   2. Constant values (either scalar or vector)
    %   3. Inline functions
    % All three may be associated with activation energy (Ea).

    numberRegex = '(?:[+\-]?(?:0|[1-9]\d*)(?:\.\d+)?(?:[eE][+\-]?\d+)?|NaN|Inf)';
    vectorRegex = sprintf('\\[(?:%s[;,])+%s[;,]?\\]',numberRegex,numberRegex);
    fcnHeadRegex = '^@\((?:(\w+),)+(\w+)\)';
    actEngeryClauseRegex = sprintf('\\.\\*exp\\((?:%s|%s)\\*\\(1/298\\.15-1/T\\)/8\\.31446\\)$',numberRegex,vectorRegex);
    actEngeryRegex = sprintf('\\.\\*exp\\((%s|%s)\\*\\(1/298\\.15-1/T\\)/8\\.31446\\)',numberRegex,vectorRegex);

    fcnStr = func2str(fcn);  % Convert function declaration to string.

    % Separate function header from body.
    idxHeadEnd = regexp(fcnStr,fcnHeadRegex,'end');
    fcnHead = fcnStr(1:idxHeadEnd);
    fcnBodyStr = fcnStr(idxHeadEnd+1:end);

    % Check for activation energy clause, separate from body if present.
    idxActEnergy = regexp(fcnBodyStr,actEngeryClauseRegex,'start');
    if isempty(idxActEnergy)
        fcnBase = fcn;
        Ea = 0;
    else
        fcnEaStr = fcnBodyStr(idxActEnergy:end);
        fcnBodyStr = fcnBodyStr(1:idxActEnergy-1);
        fcnBase = eval([fcnHead fcnBodyStr]);
        % Extract activation energy from string.
        Ea = regexp(fcnEaStr,actEngeryRegex,'tokens');
        Ea = eval(Ea{1}{1}); % Could be a vector or scalar at this point.
    end

    % Evalulate function to determine output size.
    output = fcnBase(0,0);

    % Determine if the function implements a LUT. If so, extract the
    % LUT vectors.
    isLUTParameter = false;
    if startsWith(fcnBodyStr,'interp1') || startsWith(fcnBodyStr,'(interp1')
        vectorsLUT = regexp(fcnBodyStr,vectorRegex,'match');
        if length(vectorsLUT) < 2
            warning('Failed to parse LUT from: %s',fcnStr);
        else
            vectorsLUT = cellfun(@eval,vectorsLUT,'UniformOutput',false);
            isLUTParameter = true;
        end
    end

    % Extract function parameter names from header, determine if any appear
    % in the function body. If any do, then this parameter was specified as
    % an inline function.
    isFcnParameter = false;
    fcnParams = regexp(fcnHead,'\w+','match');
    for k = 1:length(fcnParams)
        param = fcnParams{k};
        paramRegex = sprintf('\\W+%s\\W+',param);
        match = regexp(fcnBodyStr,paramRegex,'once');
        if ~isempty(match)
            isFcnParameter = true;
            break;
        end
    end

    % Structure the output data.
    if isLUTParameter
        obj.type = 'LUT';
        % Third vector may be temperature coefficients for MSMR model, not
        % needed; use only the first two.
        obj.value = vectorsLUT(1:2);
        % Use LUT value at 25degC.
        obj.value{2} = fcn(obj.value{1},273.15+25);
    elseif isFcnParameter
        obj.type = 'function';
        obj.value = fcnBase;
    elseif isscalar(output)
        obj.type = 'scalar';
        obj.value = output;
    else
        obj.type = 'vector';
        obj.value = output;
    end
    obj.Ea = Ea;
end