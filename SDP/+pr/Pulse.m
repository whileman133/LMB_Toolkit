classdef Pulse < handle
    %PULSE Data collected during a single pulse in the pulse-resistance study.

    properties(SetAccess=protected)
        time       % Time vector [s]
        voltage    % Voltage vector [V]
        current    % Current vector [A]
        soc        % Cell state-of-charge to which pulse corresponds (what was aimed for, not calculated - see QdisAh) [%]
        ipulse     % Pulse magnitude to which pulse corresponds [A]
        index      % Numerical index indicating the pulse repetition, 1-based
        QdisAh     % Total amount of charge removed from cell before start of the pulse [Ah].

        t0         % Time at pulse onset [s].
        tf         % Time at pulse offset [s].
        fs         % Nominal sampling rate [Hz].

        % Flag used to indicate empty pulse objects (needed when storing
        % pulses in arrays with uniform size).
        nil
    end

    methods
        function obj = Pulse(time, voltage, current, soc, ipulse, index, QdisAh, t0, tf, fs)
            if nargin == 0
                obj.nil = true;
            else
                obj.nil = false;
                obj.time = time;
                obj.voltage = voltage;
                obj.current = current;
                obj.soc = soc;
                obj.ipulse = ipulse;
                obj.index = index;
                obj.QdisAh = QdisAh;
                obj.t0 = t0;
                obj.tf = tf;
                obj.fs = fs;
            end
        end
    end
end