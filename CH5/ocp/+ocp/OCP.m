classdef OCP < handle
    %OCP Common interface for OCP curves.

    properties(SetAccess=protected)
        vmin      % OCP at which SOC=0% [V]
        vmax      % OCP at which SOC=100% [V]
        temp      % Temperature at which OCP data was taken [degC]
        ZU        % Relative stiochimetry vector (for potential vector) [unitless]
        U         % Potential [V]
        ZdU       % Relative stiochimetry vector (for differential potential vector) [unitless]
        dU        % Differential potential, dU/dz [V]
        meta      % Metadata string, indicates how OCP curve was derived

        % Computed from the above.
        zmin      % Relative stiochiometry where OCP=vmax [unitless]
        zmax      % Relative stiochiometry where OCP=vmin [unitless]
    end

    methods
        function obj = OCP(vmin, vmax, temp, ZU, U, ZdU, dU, meta)
            obj.vmin = vmin;
            obj.vmax = vmax;
            obj.temp = temp;
            obj.ZU = ZU;
            obj.U = U;
            obj.ZdU = ZdU;
            obj.dU = dU;
            if exist('meta','var')
                obj.meta = meta;
            end

            % Compute zmin/zmax.
            obj.zmin = obj.invUocp(obj.vmax);
            obj.zmax = obj.invUocp(obj.vmin);
        end

        function U = Uocp(obj, z)
            U = interp1(obj.ZU,obj.U,z,'linear','extrap');
        end

        function dU = dUocp(obj, z)
            dU = interp1(obj.ZdU,obj.dU,z,'linear','extrap');
        end

        function z = invUocp(obj, U)
            z = fzero(@(x)(obj.Uocp(x)-U),[0 1]);
        end
    end
end