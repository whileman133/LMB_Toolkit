function params = getFixedParams(modelspec,varargin)
    params = fastopt.unpack([],modelspec,'fixedOnly',true,varargin{:});
end