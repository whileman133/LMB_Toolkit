classdef Script < handle
    %OCPSCRIPT Stores data collected during an OCP script.

    properties(SetAccess=protected)
        name            % Name of the script
        nbr             % Number / index of the script
        temperature     % Temperature at which script was performed [C]

        time            % Time vector [s]
        current         % Current vector, signed (positive = discharge) [A]
        voltage         % Voltage vector [V]
        step            % Vector, indicates substep in which span of data was collected

        disCurrent      % Discharge current vector, always >= 0 [A]
        chgCurrent      % Charge current vector, always >= 0 [A]

        chgAh           % Cumulative charge put into cell since start of script [Ah]
        disAh           % Cumulative charge removed from cell since start of script [Ah]
    end

    methods
        function obj = Script(name, number, temperature, time, current, voltage, step)
            %OCPSCRIPT Construct an OCP script from time, current, voltage,
            %  and step vectors.

            obj.name = name;
            obj.nbr = number;
            obj.temperature = temperature;
            obj.time = time;
            obj.current = current;
            obj.voltage = voltage;
            obj.step = step;

            % Separate out discharge and charge current vectors.
            dis = current >= 0;                  chg = ~dis;
            disCurrent = zeros(size(current));   chgCurrent = zeros(size(current));
            disCurrent(dis) = abs(current(dis)); chgCurrent(chg) = abs(current(chg));
            obj.disCurrent = disCurrent;
            obj.chgCurrent = chgCurrent;

            % Calculate cumulative charge using trapaziodal method.
            obj.disAh = cumtrapz(obj.time, disCurrent)/3600;
            obj.chgAh = cumtrapz(obj.time, chgCurrent)/3600;
        end

        function s = toStruct(obj)
            s.name = obj.name;
            s.nbr = obj.nbr;
            s.temperature = obj.temperature;
            s.time = obj.time;
            s.current = obj.current;
            s.voltage = obj.voltage;
            s.step = obj.step;
            s.disCurrent = obj.disCurrent;
            s.chgCurrent = obj.chgCurrent;
            s.chgAh = obj.chgAh;
            s.disAh = obj.disAh;
        end
    end
end