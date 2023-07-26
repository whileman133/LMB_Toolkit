function model = unpack(vect, metadata, varargin)
    %UNPACK Unstuff a model structure from a vector.
    %
    % model = UNPACK(vect, metadata) converts the vector VECT back into
    %   a model structure using the metadata METADATA.
    %
    % model = UNPACK(...,'sparse',true) omits fixed parameter values from
    %   the model.
    %
    % model = UNPACK(...,'flat',true) returns a flat-parameter model
    %   instead of a hierarchical model.

    parser = inputParser;
    parser.addOptional('sparse',false,@islogical);
    parser.addOptional('flat',false,@islogical);
    parser.parse(varargin{:});
    sparse = parser.Results.sparse;
    flat = parser.Results.flat;

    flatmodel = struct;
    vect = vect(:);

    % Extract parameters.
    paramnames = fieldnames(metadata.params);
    cursor = 1;
    for k = 1:length(paramnames)
        paramname = paramnames{k};
        meta = metadata.params.(paramname);

        % Determine multiplicity.
        mult = 1;
        if strcmpi(meta.tempfcn,'lut')
           mult = metadata.ntemps;
        end

        % Determine length.
        len = meta.len;
        if isfield(meta,'fix')
            len = sum(meta.fixmask);
        end

        % Extract value and store into model structure.
        value = vect(cursor:cursor+len*mult-1);
        value = reshape(value,[len,mult]);
        if isfield(meta,'fix')
            if ~sparse
                flatmodel.(paramname) = meta.fix;
                flatmodel.(paramname)(meta.fixmask) = value;
            elseif len > 0
                flatmodel.(paramname) = value;
            else
                % Fixed parameter; omit from sparse output.
            end
        else
            flatmodel.(paramname) = value;
        end % if
        cursor = cursor + len*mult;

        % Decode logarithmic values.
        if isfield(meta,'logscale') && meta.logscale == true && isfield(flatmodel,paramname)
            flatmodel.(paramname) = 10.^(flatmodel.(paramname));
        end

        % Check for activation energy.
        if strcmpi(meta.tempfcn,'Eact')
            % Activation energy is present; extract from vector.
            paramnameEact = [paramname 'Eact'];
            Eact = vect(cursor);
            flatmodel.(paramnameEact) = Eact;
            cursor = cursor + 1;
        end % if
    end % for

    if flat
        model = flatmodel;
    else
        % Deflatten model.
        model = struct;
        paramnames = fieldnames(flatmodel);
        for k = 1:length(paramnames)
            paramname = paramnames{k};
            parts = split(paramname,'__');
            model = setfield(model,parts{:},flatmodel.(paramname));
        end
    end
end