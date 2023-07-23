classdef Test < handle
    %TEST Data collected for a pulse-resistance test.

    properties(SetAccess=protected)
        name     % Name of the test (cell name for example.)
        temp     % Temperature at which test was performed [C]
        pulses   % 3-D structure array of Pulse objects.
                 %   1st dimension: SOC
                 %   2nd dimension: pulse current
                 %   3rd dimension: pulse repetitions
                 % Note: the 3rd dimension is ragged; "empty" pulse objects
                 % are indicated by a true value of the `nil` property.
                 % Use `getPules()` to simplify fetching nonempty pulses.

        repeat   % Maximum number of repetitions for each pulse.
        socs     % SOC values at which the test was performed (what was aimed for, not calculated) [%]
        currents % Current magnitudes at which the test was performed [A]
        QdisAh   % Amount of charge removed from the cell at each SOC setpoint [Ah].
                 % Use to calculate the true SOC of the cell: DOD = QdisAh./Q
                 % where Q is the capacity of the cell, SOC = 1 - DOD.
    end

    methods
        function obj = Test(name, temp, pulseVect)
            %PRTEST Create a pulse-resistance test from a vector of Pulse 
            %  objects.

            obj.name = name;
            obj.temp = temp;

            % First determine the SOC and current vectors, and the maximum
            % number of pulse repetitions.
            pulseVect = pulseVect(:);
            repeat = max([pulseVect.index]);
            socs = unique([pulseVect.soc]);         % Note that unique.m also sorts ascending!
            currents = unique([pulseVect.ipulse]);
            QdisAh = unique([pulseVect.QdisAh]);
            QdisAh = fliplr(QdisAh);                % DOD needs to be descending to be consistent with ascending SOC.

            % Store each pulse in a 3-D structure array.
            pulses = pr.Pulse.empty(length(socs), length(currents), 0);
            for idxPulse = 1:length(pulseVect)
                pulse = pulseVect(idxPulse);
                [~, idxSOC] = min(abs(socs - pulse.soc));
                [~, idxCurrent] = min(abs(currents - pulse.ipulse));
                idxRepeat = pulse.index;
                pulses(idxSOC,idxCurrent,idxRepeat) = pulse;
            end

            % Store results to object.
            obj.repeat = repeat;
            obj.socs = socs;
            obj.currents = currents;
            obj.QdisAh = QdisAh;
            obj.pulses = pulses;
        end

        function pulses = getPulses(obj, soc, current, tol)
            %GETPULSES Fetch all pulses at the specified SOC and current.

            if ~exist('tol','var')
                tol = 1e-9;
            end

            [errSOC, idxSOC] = min(abs(obj.socs - soc));
            [errCurrent, idxCurrent] = min(abs(obj.currents - current));

            if errSOC > tol
                error('Matching SOC not found.');
            end
            if errCurrent > tol
                error('Matching current not found.');
            end
            
            pulses = obj.getPulsesByIndex(idxSOC, idxCurrent);
        end

        function pulses = getPulsesByIndex(obj, idxSOC, idxCurrent)
            %GETPULSESBYINDEX Fetch all pulses at the specified SOC and
            %  current indicies.

            pulses = obj.pulses(idxSOC,idxCurrent,:);
            pulses = pulses(:);                     % Convert 3D structure to column vector.
            pulses = pulses([pulses.nil] == false); % Remove empty pulses.
            pulses = pulses';                       % Convert to row vector.
        end

        function s = shallow(obj)
            %SHALLOW Return a shallow copy of this test without the
            % underlying pulse data but with everything else.

            fields = {'name','temp','repeat','socs','currents','QdisAh'};
            for k = 1:length(fields)
                field = fields{k};
                s.(field) = obj.(field);
            end
        end
    end
end