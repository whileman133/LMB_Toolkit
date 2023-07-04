function model = unpack(vect, metadata, varargin)
    %UNPACK Unstuff a model from a vector.
    %
    % model = UNPACK(vect, metadata) converts the vector VECT back into
    %   a model using the metadata METADATA.
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
        if isfield(meta,'fix')
            len = sum(meta.fixmask);
            if ~sparse
                flatmodel.(paramname) = meta.fix;
                flatmodel.(paramname)(meta.fixmask) = vect(cursor:cursor+len-1);
            elseif len > 0
                flatmodel.(paramname) = vect(cursor:cursor+len-1);
            end
            cursor = cursor + len;
        else
            flatmodel.(paramname) = vect(cursor:cursor+meta.len-1);
            cursor = cursor + meta.len;
        end
        if isfield(meta,'logscale') && meta.logscale == true && isfield(flatmodel,paramname)
            % Decode logarithmic value.
            flatmodel.(paramname) = 10.^(flatmodel.(paramname));
        end
    end

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