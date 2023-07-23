classdef Surface < handle
    %PREstimate Pulse-resistance surface estimate derived from
    % laboratory measurements.

    properties
        name        % Name of the test from which estimate was derived (cell name for example.)
        temp        % Temperature of test from which estimate was derived [C]
        QAh         % Total capacity of the cell [Ah]
        modelType   % Type of ECM used (class name)
        opts        % Structure of options used in fitting ECM to laboratory data
        z           % Cell state-of-charge vector (true SOC, not what was programmed) [fractional]
        socs        % Cell state-of-charge vector (what was programmed, not exact) [%]
        iapp        % Pulse current magnitude vector [A]
        models      % 2D cell array of fit ECM arrays: first index <=> SOC, second index <=> iapp
                    % (multiple ECMs per SOC/iapp setpoint due to pulse repetitions)
    end

    methods
        function obj = Surface(prtest, ocptest, models, opts)
            obj.name = prtest.name;
            obj.temp = prtest.temp;
            obj.QAh = ocptest.QAh;
            obj.z = 1 - prtest.QdisAh/ocptest.QAh;
            obj.socs = prtest.socs;
            obj.iapp = prtest.currents;
            obj.models = models;
            obj.opts = opts;
            obj.modelType = class(models{1,1});
        end

        function [P, Pbound] = getParam(obj, ECMfcn, varargin)
            %GETPARAM Fetch the specified ECM parameter at each SOC/iapp
            % setpoint.

            if ~isa(ECMfcn,'function_handle')
                % Assume ECMfcn specifies the string name of an ECM parameter.
                ECMfcn = @(model)model.(ECMfcn);
            end

            parser = inputParser;
            parser.addOptional('removeOutliers',false,@islogical);
            parser.parse(varargin{:});
            rmoutliers = parser.Results.removeOutliers;

            P = zeros(length(obj.iapp),length(obj.z));
            Pbound = zeros(size(P));
            for idxI = 1:length(obj.iapp)
                for idxZ = 1:length(obj.z)
                    ECMs = obj.models{idxI,idxZ};
                    p = arrayfun(ECMfcn,ECMs);
                    if rmoutliers
                        outliers = isoutlier(p);
                        p = p(~outliers);
                    end
                    P(idxI,idxZ) = mean(p);
                    if length(ECMs) == 1
                        try
                            Pbound(idxI,idxZ) = 3*ECMfcn(ECMs(1).bound);
                        catch
                            Pbound(idxI,idxZ) = NaN;
                        end
                    else
                        Pbound(idxI,idxZ) = 3*std(p);
                    end
                end
            end
        end % getParam()

        function ecms = getECMs(obj, soc, current, tol)
            %GETECMS Fetch all ECMs at the specified SOC and current.

            if ~exist('tol','var')
                tol = 1e-9;
            end

            [errSOC, idxSOC] = min(abs(obj.z*100 - soc));
            [errCurrent, idxCurrent] = min(abs(obj.iapp - current));

%             if errSOC > tol
%                 error('Matching SOC not found.');
%             end
%             if errCurrent > tol
%                 error('Matching current not found.');
%             end
            
            ecms = obj.models{idxCurrent,idxSOC};
        end
    end
end