function modelVect = splittemps(model,modelspec)
%SPLITTEMPS Create model structures with multiplicty=1 at each temperature.

temps = modelspec.temps;
ntemps = modelspec.ntemps;

% Flatten model.
model = fastopt.flattenstruct(model);

% Pre-allocate flat models.
flatmodels = [];
for m = ntemps:-1:1
    flatmodels(m) = model;
end

% Compute parameter values at each temperature.
paramnames = fieldnames(metadata.params);
for k = 1:length(paramnames)
    paramname = paramnames{k};
    meta = metadata.params.(paramname);
    value = model.(paramname);

    if strcmpi(meta.tempfcn,'fix')
        for m = ntemps:-1:1
            flatmodels(m).
        end
    elseif strcmpi(meta.tempfcn,'lut')
    elseif strcmpi(meta.tempfcn,'Eact')
    end
end % for

end