function [vect, metadata] = pack(model, metadata, varargin)
    %PACK Stuff a model into a vector.
    %
    % vect = PACK(model,metadata) converts the pulse model MODEL into a
    %   vector VECT for use with optimization routines. META is generated
    %   using MODELSPEC().

    parser = inputParser;
    parser.addOptional('default',[]);
    parser.addOptional('coerce',false,@islogical);
    parser.addOptional('dim',1,@isscalar);
    parser.parse(varargin{:});
    defaultval = parser.Results.default;
    coerce = parser.Results.coerce;
    dim = parser.Results.dim;

    % Flatten model.
    model = fastopt.flattenstruct(model);

    % Construct parameter vector.
    vect = zeros(metadata.nvars,dim);
    paramnames = fieldnames(metadata.params);
    cursor = 1;
    for k = 1:length(paramnames)
        paramname = paramnames{k};
        meta = metadata.params.(paramname);

        if isfield(meta,'fix') && ~any(meta.fixmask)
            % Fixed parameter, skip packing.
            continue;
        end

        % Regular parameter.
        if isfield(model,paramname)
            value = model.(paramname);
        elseif ~isempty(defaultval)
            value = repmat(defaultval,meta.len,1);
        else
            error('Could not resolve value for parameter %s',paramname);
        end
        if isfield(meta,'logscale') && meta.logscale == true
            % Encode logarithmic value.
            value = log10(value);
        end

        % Expand scalar parameters to correct size.
        [vallen, ~] = size(value);
        if coerce && vallen == 1
            value = value*ones(meta.len,dim);
        end

        % Remove fixed components of the value (if it is a vector).
        if isfield(meta,'fix')
            value = value(meta.fixmask);
            len = sum(meta.fixmask);
        else
            len = meta.len;
        end

        % Place value into parameter vector.
        [vallen, valdim] = size(value);
        if vallen == len && valdim == dim
            vect(cursor:cursor+len-1,:) = value(:,:);
        else
            error(['Parameter %s has incorrect length ' ...
                '(metadata says it should be %d, but it is actualy %d)'], ...
                paramname,len,length(value));
        end
        cursor = cursor + len;
    end % for
end