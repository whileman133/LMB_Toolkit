function hierarchy = unflattenstruct(flat)
    %UNFLATTENSTRUCT Deflatten a structure flattened by flattenstruct.
    hierarchy = struct;
    paramnames = fieldnames(flat);
    for k = 1:length(paramnames)
        paramname = paramnames{k};
        parts = split(paramname,'__');
        hierarchy = setfield(hierarchy,parts{:},flat.(paramname));
    end
end

