function [low, high] = checkbounds(modelspec,values,varargin)
%CHECKBOUNDS

p = inputParser;
p.addParameter('lb',struct,@isstruct);
p.addParameter('ub',struct,@isstruct);
p.addParameter('tol',0.02,@isscalar);
p.parse(varargin{:});

values = fastopt.flattenstruct(values);
lb = fastopt.flattenstruct(p.Results.lb);
ub = fastopt.flattenstruct(p.Results.ub);
tol = p.Results.tol;
low = struct;
high = struct;

paramnames = fieldnames(modelspec.params);
for k = 1:length(paramnames)
    name = paramnames{k};
    meta = modelspec.params.(name);
    if isfield(meta,'fix') && ~any(meta.fixmask)
        % Fixed parameter.
        continue;
    end
    value = values.(name);
    if isfield(lb,name)
        lower = lb.(name);
        l = (value-lower)./lower <= tol;
        if isfield(meta,'fix')
            l = l(meta.fixmask);
        end
        if any(l)
            low.(name) = l;
        end
    end
    if isfield(ub,name)
        upper = ub.(name);
        h = (upper-value)./upper <= tol;
        if isfield(meta,'fix')
            h = h(meta.fixmask);
        end
        if any(h)
            high.(name) = h;
        end
    end
end


end

