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
    %
    % model = UNPACK(...,'fixed',true) returns only the values of fixed
    %   parameters. In this case, the input VECT is optional and may be
    %   replaced by the empty vector [].

    parser = inputParser;
    parser.addRequired('vect',@(x)isnumeric(x)&&isvector(x)||isempty(x));
    parser.addRequired('modelspec',@(x)isstruct(x)&&isscalar(x));
    parser.addOptional('sparse',false,@islogical);
    parser.addOptional('flat',false,@islogical);
    parser.addOptional('fixedOnly',false,@islogical);
    parser.parse(vect,metadata,varargin{:});
    sparse = parser.Results.sparse;
    flat = parser.Results.flat;
    fixedOnly = parser.Results.fixedOnly;

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
        if isempty(vect)
            value = [];
        else
            value = vect(cursor:cursor+len*mult-1);
            value = reshape(value,[len,mult]);
        end
        if isfield(meta,'fix')
            if ~sparse
                flatmodel.(paramname) = meta.fix;
                flatmodel.(paramname)(meta.fixmask) = value;
            elseif len > 0
                % Some components are not fixed; omit the fixed components
                % and retain the unfixed components.
                flatmodel.(paramname) = value;
            else
                % Fixed parameter; omit from sparse output.
            end
        elseif ~fixedOnly
            % Unfixed parameter.
            flatmodel.(paramname) = value;
        end % if
        cursor = cursor + len*mult;

        % Decode logarithmic values.
        if isfield(meta,'logscale') && meta.logscale == true && isfield(flatmodel,paramname)
            flatmodel.(paramname) = 10.^(flatmodel.(paramname));
        end

        % Check for activation energy.
        if strcmpi(meta.tempfcn,'Eact') && ~isempty(vect)
            % Activation energy is present; extract from vector.
            paramnameEact = [paramname '_Eact'];
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