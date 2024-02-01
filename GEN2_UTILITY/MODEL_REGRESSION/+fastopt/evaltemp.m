function modelVect = evaltemp(model,modelspec,TdegC)
%EVALTEMP Create model structures of multiplicty=1 at each temperature.

temps = TdegC+273.15;
ntemps = length(TdegC);
Tref = modelspec.Tref;
R = TB.const.R;

% Flatten model.
flatmodel = fastopt.flattenstruct(model);

% Pre-allocate model vector.
clear flatmodels;
for m = ntemps:-1:1
    flatmodels(m) = flatmodel;
end

% Compute parameter values at each temperature.
paramnames = fieldnames(modelspec.params);
for k = 1:length(paramnames)
    paramname = paramnames{k};
    meta = modelspec.params.(paramname);
    if isfield(meta,'fix') && ~any(meta.fixmask) && ~isfield(flatmodel,paramname)
        % OK to omit fixed parameters
        continue;
    end
    value = flatmodel.(paramname);
    
    if strcmpi(meta.tempfcn,'fix')
        % Same value for all temperatures.
        for m = ntemps:-1:1
            flatmodels(m).(paramname) = value;
        end
    elseif strcmpi(meta.tempfcn,'lut')
        % Lookup table vs temperature.
        for m = ntemps:-1:1
            flatmodels(m).(paramname) = value(:,m);
        end
    elseif strcmpi(meta.tempfcn,'Eact')
        % Ahrrenius relation vs temperature. 
        % Assumes value is reference value at Tref.
        paramnameEact = [paramname '_Eact'];
        Eact = flatmodel.(paramnameEact);
        for m = ntemps:-1:1
            T = temps(m);
            if strcmpi(meta.tempcoeff,'+')
                flatmodels(m).(paramname) = value.*exp((Eact/R)*(1/Tref-1/T));
            else
                flatmodels(m).(paramname) = value./exp((Eact/R)*(1/Tref-1/T));
            end
        end
        % Remove Eact.
        flatmodels = rmfield(flatmodels,paramnameEact);
    else
        error('Unrecognized temperature function: %s',meta.tempfcn);
    end % if
end % for

% Unflatten models.
clear modelVect;
for m = ntemps:-1:1
    modelVect(m) = fastopt.unflattenstruct(flatmodels(m));
end

end