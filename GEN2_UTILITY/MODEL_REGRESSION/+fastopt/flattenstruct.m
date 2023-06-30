function flat = flattenstruct(s, prefix, whitelist)
    %FLATTENSTRUCT Recursively flatten a structure with hierarchy.
    %
    % flat = FLATTENSTRUCT(s) converts the hierarchical structure S
    %   into flat structure FLAT by joining nested field names with
    %   the double underscore (__) delimiter.

    flat = struct;

    if ~exist('prefix','var')
        prefix = '';
    end
    if ~exist('whitelist','var')
        whitelist = {};
        hasWhitelist = false;
    else
        hasWhitelist = true;
    end
    if ~isstruct(s)
        flat.(prefix) = s;
        return;
    end

    fnames = fieldnames(s);
    if any(strcmp('noflatten__',fnames))
        flat.(prefix) = s;
        return;
    end

    for k = 1:length(fnames)
        fname = fnames{k};
        if hasWhitelist && ~any(strcmp(fname,whitelist))
            % Skip this field.
            continue;
        end
        value = s.(fname);
        if strcmp(prefix,'')
            flatName = fname;
        else
            flatName = [prefix '__' fname];
        end
        sflat = fastopt.flattenstruct(value,flatName);
        subfnames = fieldnames(sflat);
        for j = 1:length(subfnames)
            subfname = subfnames{j};
            flat.(subfname) = sflat.(subfname);
        end
    end % for
end