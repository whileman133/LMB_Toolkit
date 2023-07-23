classdef PulseModel < handle
    %PULSEMODEL Parameters for a LMB pulse model.

    properties(SetAccess=protected)
        neg   % Negative electrode interface model.
        pos   % Positive electrode interface model.

        % Region conductances (lumped model).
        kd   % Electrolyte conductance, dead lithium layer.
        ks   % Electrolyte conductance, separator.
        kp   % Electrolyte conductance, positive electrode.
        sp   % Solid conductance, positive electrode.
    end

    methods
        function obj = PulseModel(neg, pos, kd, ks, kp, sp)
            %PULSEMODEL
            obj.neg = neg;
            obj.pos = pos;
            obj.kd = kd;
            obj.ks = ks;
            obj.kp = kp;
            obj.sp = sp;
        end
    end
end