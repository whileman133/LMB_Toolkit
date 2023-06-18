function m = modelspec(params)
    %MODELSPEC Generate a description of a model.
    %
    % m = MODELSPEC(params) creates metadata structure M for 
    %   packing and unpacking a model. PARAMS is a structure 
    %   describing the parameters that make up the model. Use PARAM() to 
    %   generate values to assign to fields of the PARAMS structure.
    %   Fields of PARAMS may be nested arbitrarily.

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
            nvars = nvars + len;
        else
            % Packed variable; increment degrees-of-freedom counter.
            len = paramspec.len;
            nvars = nvars + len;
        end
        m.params.(paramname) = paramspec;
    end
    m.nvars = nvars;
end