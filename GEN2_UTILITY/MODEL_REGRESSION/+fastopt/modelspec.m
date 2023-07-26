function m = modelspec(params,varargin)
    %MODELSPEC Generate a description of a model.
    %
    % m = MODELSPEC(params) creates metadata structure M for 
    %   packing and unpacking a model. PARAMS is a structure 
    %   describing the parameters that make up the model. Use PARAM() to 
    %   generate values to assign to fields of the PARAMS structure.
    %   Fields of PARAMS may be nested arbitrarily, as the structure is
    %   flattened by this function.

    parser = inputParser;
    parser.addRequired('params',@isstruct);
    parser.addParameter('tempsdegC',25,@(x)isnumeric(x)&&isvector(x)&&length(x)>=1);
    parser.addParameter('TrefdegC',25,@(x)isscalar(x));
    parser.parse(params,varargin{:});
    arg = parser.Results;
    ntemps = length(arg.tempsdegC);

    % Flatten parameter structure.
    params = fastopt.flattenstruct(params);

    % Build metadata.
    paramnames = fieldnames(params);
    nparams = length(paramnames);
    nvars = 0;
    for k = 1:nparams
        paramname = paramnames{k};
        paramspec = params.(paramname);
        if isfield(paramspec,'fix')
            % Fixed variable (not packed into vector); do not increment 
            % degrees-of-freedom counter.
            % Exception: if this is a vector and some of the components are
            % not fixed, then they get packed into the vector and the
            % degrees-of-freedom counter is incremented.
            len = sum(paramspec.fixmask);
        else
            % Packed variable; increment degrees-of-freedom counter.
            len = paramspec.len;
        end % if
        if len>0
            if strcmpi(paramspec.tempfcn,'lut')
                % Add additional parameters for each temperature.
                len = len*ntemps;
            elseif strcmpi(paramspec.tempfcn,'Eact')
                % Add additional parameter for activation energy.
                len = len + 1;
            end % if
        end % if
        nvars = nvars + len;
        m.params.(paramname) = paramspec;
    end
    m.nvars = nvars;
    m.temps = arg.tempsdegC+273.15;
    m.ntemps = ntemps;
    m.Tref = arg.TrefdegC+273.15;
end