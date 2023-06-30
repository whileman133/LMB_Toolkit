classdef MSMRFit < handle
    %MSMRFIT Regressed MSMR model derived from laboratory measurements.

    properties(SetAccess=protected)
        name        % Name of the test from which fit model was derived (cell name for example.)
        temp        % Temperature of test from which fit model was derived [C]
        ocpest      % OCP estimate from which the fit model was derived
        opts        % Structure of options used in fitting ECM to laboratory data
        model       % Regressed MSMR model (instance of MSMR).
    end

    methods
        function obj = MSMRFit(ocptest, ocpest, model, opts)
            obj.name = ocptest.name;
            obj.temp = ocptest.temp;
            obj.ocpest = ocpest;
            obj.model = model;
            obj.opts = opts;
        end

        function s = toStruct(obj)
            s.name = obj.name;
            s.temp = obj.temp;
            s.opts = obj.opts;
            s.MSMR.U0 = obj.model.Uj0;
            s.MSMR.X = obj.model.Xj;
            s.MSMR.omega = obj.model.Wj;
            s.MSMR.zmin = obj.model.zmin;
            s.MSMR.zmax = obj.model.zmax;
            s.MSMR.J = length(obj.model.Xj);
            s.ocptest = obj.ocpest.ocptest.toStruct();
        end
    end
end