classdef ECMFit < handle
    %ECMFit Data associated with an Equivalent Circuit Model (ECM)
    %  regression to pulse-response measurements.

    properties(SetAccess=protected)
        model     % Regressed ECM
        opts      % Structure of ECM-specific options used in the regression
        soc       % Cell state-of-charge [%]
        iapp      % Pulse current magnitude [A]
    end

    methods
        function obj = ECMFit(varargin)
            p = inputParser;
            addOptional(p,'model',nan);
            addOptional(p,'opts',nan,@isstruct);
            addOptional(p,'soc',nan);
            addOptional(p,'iapp',nan);
            parse(p,varargin{:});

            params = struct2cell(p.Results);
            idxnan = cellfun(...
                @(x)(~isa(x,'function_handle') && any(isnan(x))),params);
            if any(idxnan)
                paramnames = fieldnames(p.Results);
                error('Missing required parameters %s',...
                      strjoin(paramnames(idxnan),', '));
            end

            obj.model = p.Results.model;
            obj.opts = p.Results.opts;
            obj.soc = p.Results.soc;
            obj.iapp = p.Results.iapp;
        end
    end
end