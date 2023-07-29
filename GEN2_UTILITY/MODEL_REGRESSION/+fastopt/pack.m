function [vect, metadata] = pack(model,metadata,varargin)
    %PACK Stuff a model structure into a vector.
    %
    % vect = PACK(model,metadata) converts the model structure MODEL into
    %   vector VECT for use with optimization routines. META is generated
    %   using MODELSPEC().

    parser = inputParser;
    parser.addOptional('default',[]);
    parser.addOptional('coerce',true,@islogical);
    parser.parse(varargin{:});
    defaultval = parser.Results.default;
    coerce = parser.Results.coerce;

    % Flatten model.
    model = fastopt.flattenstruct(model);

    % Construct parameter vector.
    vect = zeros(metadata.nvars,1);
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
        [vallen, valmult] = size(value);

        % Determine required multiplicity of the parameter (>1 when 
        % temperature dependence is modeled using LUT approach).
        mult = 1;
        if strcmpi(meta.tempfcn,'lut')
           mult = metadata.ntemps;
        end

        % Expand scalar parameters to correct size.
        if coerce && vallen == 1 && valmult == 1
            value = value*ones(meta.len,mult);
        end

        % Remove fixed components of the value (if it is a vector).
        if isfield(meta,'fix')
            value = value(meta.fixmask);
            len = sum(meta.fixmask);
        else
            len = meta.len;
        end

        % Place value into parameter vector.
        [vallen, valmult] = size(value);
        if vallen == len && valmult == mult
            vect(cursor:cursor+len*mult-1,1) = value(:);
        else
            if vallen ~= len
                error(['Parameter %s has incorrect length ' ...
                    '(metadata says it should be %d, but it is actualy %d)'], ...
                    paramname,len,vallen);
            else
                error(['Parameter %s has incorrect multiplicity ' ...
                    '(metadata says it should be %d, but it is actualy %d)'], ...
                    paramname,mult,valmult);
            end
        end
        cursor = cursor + len*mult;

        % Check for activation energy parameter.
        if strcmpi(meta.tempfcn,'Eact')
            % Should exist - add Eact to parameter vector.
            paramnameEact = [paramname '_Eact'];
            if ~isfield(model,paramnameEact)
                error(['Eact not found for parameter %s. ' ...
                    'Hint: Specify an activation energy parameter %s.'], ...
                    paramname,paramnameEact);
            end
            Eact = model.(paramnameEact);
            [Eactlen, Eactmult] = size(Eact);
            if Eactlen ~= 1
                error('%s must have length 1.',paramnameEact);
            end
            if Eactmult ~= 1
                error('%s must have multiplicity 1.',paramnameEact);
            end
            vect(cursor:cursor+1,1) = Eact(:);
            cursor = cursor + 1;
        end
    end % for
end