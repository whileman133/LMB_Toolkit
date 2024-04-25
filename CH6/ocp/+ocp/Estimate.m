classdef Estimate < handle
    %OCPESTIMATE Data related to an OCP estimate.

    properties(SetAccess=protected)
        ocptest   % OCP test data from which estimate was derived.
        dvbin     % Bin width for differential capacity computation.
        Z         % Relative stiochimetry [unitless]
        V         % Potential [V]
        dzdvRefZ  % Differential capacity (over stiochiometry) [1/V]
        refZ      % Stiochiometry vector for differential capacity [unitless].
        dzdvRefV  % Differential capacity (over potential) [1/V]
        refV      % Potential vector for differential capacity [V].
        meta      % Metadata string.
    end

    methods
        function obj = Estimate(ocptest, dvbin, Z, V, meta)
            obj.ocptest = ocptest;
            obj.dvbin = dvbin;
            obj.Z = Z;
            obj.V = V;
            if exist('meta','var')
                obj.meta = meta;
            end

            % Compute differential capacity.
            [obj.dzdvRefV, obj.refV, obj.dzdvRefZ, obj.refZ] = ... 
                ocp.smoothdiff(Z, V, dvbin);
        end

        function ocpObj = asOCP(obj)
            %ASOCP Convert this OCP estimate to a standard OCP object.

            ocpObj = ocp.OCP( ...
                obj.ocptest.vmin, obj.ocptest.vmax, obj.ocptest.temp, ...
                obj.Z, obj.V, ...
                obj.refZ, -abs(1./obj.dzdvRefZ), ...
                sprintf('OcpEstimate %s (OcpTest: %s)',obj.meta,obj.ocptest.name));
        end
    end
end