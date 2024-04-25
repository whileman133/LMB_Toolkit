classdef DirectEstimate < handle
    %DIRECTOCPESTIMATE Use dis/charge data to approximate OCP.

    properties(SetAccess=protected)
        ocptest   % OCP test object from which estimate is derived.
        dvbin     % Size of OCP bins used in differential capacity calculation [V].
    end

    methods
        function obj = DirectEstimate(ocptest, dvbin)
            obj.ocptest = ocptest;
            obj.dvbin = dvbin;
        end

        function est = useDis(obj)
            %USEDIS Build an OcpEstimate using discharge data.
            est = ocp.Estimate(obj.ocptest,obj.dvbin,obj.ocptest.disZ,obj.ocptest.disV,'Direct, Discharge');
        end

        function est = useChg(obj)
            %USEDIS Build an OcpEstimate using discharge data.
            est = ocp.Estimate(obj.ocptest,obj.dvbin,obj.ocptest.chgZ,obj.ocptest.chgV,'Direct, Discharge');
        end
    end
end