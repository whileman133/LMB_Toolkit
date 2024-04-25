classdef Test < handle
    %OCPTEST Data collected during an OCP test on a particular cell at a 
    %   particular temperature.

    properties(Constant)
        SCRIPTNAMES = {'script1', 'script2', 'script3', 'script4'};
    end

    properties(SetAccess=protected)
        name     % Name of the test (cell name for example.)
        temp     % Temperature at which test was performed [C]
        vmin     % Minimum cell voltage [V]
        vmax     % Maximum cell voltage [V]

        % Data collected during scripts.
        script1  % Constant-current discharge from Vmax (Temperature=T)
        script2  % Calibration step; dither, CC/CV discharge (Temperature=25C)
        script3  % Constant-current charge from Vmin (Temperature=T)
        script4  % Calibration step; dither, CC/CV charge (Temperature=25C)

        % Vector of scripts broken out above (simplifies access in some cases.)
        % Ordered according to the SCRIPTNAMES constant above.
        scripts

        % Concatenated vectors.
        % Note: at edges between scripts, we may get duplicate time values;
        %   use 'script' vector to tell them apart.
        % Note: we break out Ah charged at Tref vs Ttest since the
        %   quantities are needed in the coulombic efficiency calculation.
        time        % Time [s]
        voltage     % Cell voltage [V]
        current     % Cell current, positive for discharge [A]
        chgAhTref   % Cumulative Ah charged at reference temperature.
        chgAhTtest  % Cumulative Ah charged at test temperature.
        chgAh       % Cumulative Ah charged at either temperature.
        disAh       % Cumulative Ah discharged at either temperature.
        script      % Script number.
        temperature % Temperature at which data was collected [degC].

        % Processed data.
        eta         % Coulombic efficiency [unitless].
        QAh         % Relative (or usable) capacity [Ah].
        disZ        % Relative stiochiometry (theta) vector, discharge [unitless].
        disV        % Discharge voltage vector [V].
        chgZ        % Relative stiochiometry (theta) vector, charge [unitless].
        chgV        % Charge voltage vector [V].
    end

    methods
        function obj = Test(name, temperature, vmin, vmax, scripts)
            %OCPTEST Construct an OCP test from component scripts.

            obj.name = name;
            obj.temp = temperature;
            obj.vmin = vmin;
            obj.vmax = vmax;

            % Check scripts.
            missing = setdiff(obj.SCRIPTNAMES, {scripts.name});
            extra = setdiff({scripts.name}, obj.SCRIPTNAMES);
            if ~isempty(missing)
                string = sprintf('%s ', missing{:});
                warning('Missing scripts: %s.', string);
            end
            if ~isempty(extra)
                string = sprintf('%s ', extra{:});
                warning('Got extra (unrecognized) scripts: %s.', string);
            end

            % Assign recognized scripts to properties.
            obj.scripts = [];
            for scriptname = obj.SCRIPTNAMES
                matches = strcmp({scripts.name}, scriptname);
                if any(matches)
                    script = scripts(matches);
                    obj.(script.name) = script;
                    obj.scripts = [obj.scripts script];
                end
            end

            % Concatenate vectors.
            % Note: at edges between scripts, we may get duplicate time values.
            % Use 'script' vector to tell them apart.
            time = 0; chgAh = 0; chgAhTtest = 0; chgAhTref = 0; disAh = 0; voltage = []; current = [];
            scriptnbr = []; temperature = [];
            for script = obj.scripts
                time = [time; time(end)+script.time];
                chgAh = [chgAh; chgAh(end)+script.chgAh];
                disAh = [disAh; disAh(end)+script.disAh];
                voltage = [voltage; script.voltage];
                current = [current; script.current];
                scriptnbr = [scriptnbr; ones(size(script.time))*script.nbr];
                temperature = [temperature; ones(size(script.time))*script.temperature];

                if script.temperature == obj.temp
                    % Performed at Ttest.
                    chgAhTtest = [chgAhTtest; chgAhTtest(end)+script.chgAh];
                    chgAhTref = [chgAhTref; ones(size(script.time))*chgAhTref(end)]; % no increase in charged Ah at Tref 
                else
                    % Performed at Tref.
                    chgAhTref = [chgAhTref; chgAhTref(end)+script.chgAh];
                    chgAhTtest = [chgAhTtest; ones(size(script.time))*chgAhTtest(end)]; % no increase in charged Ah at Ttest 
                end
            end
            obj.time = time(2:end);
            obj.chgAhTref = chgAhTref(2:end);
            obj.chgAhTtest = chgAhTtest(2:end);
            obj.chgAh = chgAh(2:end);
            obj.disAh = disAh(2:end);
            obj.voltage = voltage;
            obj.current = current;
            obj.script = scriptnbr;
            obj.temperature = temperature;
        end % OcpTest()

        function processOCP(obj, eta, Q)
            %PROCESSOCP Compute dis/charge voltage vs relative
            %  stiochiometry and store the results internally.
            % 
            % Requires externally-computed values for eta and Q.

            obj.eta = eta;
            obj.QAh = Q;

            % Limit the data we consider to instants when the cell is
            % actually being dis/charged and not resting.
            idxDis = obj.script1.current ~= 0; 
            idxChg = obj.script3.current ~= 0;

            obj.disZ = (obj.script1.disAh(idxDis) - eta*obj.script1.chgAh(idxDis)) / Q;
            obj.disV = obj.script1.voltage(idxDis);
            obj.chgZ = 1 + (obj.script3.disAh(idxChg) - eta*obj.script3.chgAh(idxChg)) / Q;
            obj.chgV = obj.script3.voltage(idxChg);

            % Due to floating point arithmetic error in numerical integration, 
            % sometimes disZ does not start precisely at 0 and chgZ does not 
            % start precisely at 1; subtract off the differences.
            obj.disZ = obj.disZ - obj.disZ(1);
            obj.chgZ = obj.chgZ - obj.chgZ(1) + 1;

            % The charge vectors are backwards due to the way the test is
            % performed; flip for consistency with discharge vectors (i.e.,
            % low stiochiometry to high stiochiometry.)
            obj.chgZ = flip(obj.chgZ);
            obj.chgV = flip(obj.chgV);
        end

        function s = toStruct(obj)
            s.name = obj.name;
            s.temp = obj.temp;
            s.vmin = obj.vmin;
            s.vmax = obj.vmax;
            s.script1 = obj.script1.toStruct();
            s.script2 = obj.script2.toStruct();
            s.script3 = obj.script3.toStruct();
            s.script4 = obj.script4.toStruct();
            s.time = obj.time;
            s.voltage = obj.voltage;
            s.current = obj.current;
            s.chgAhTref = obj.chgAhTref;
            s.chgAhTtest = obj.chgAhTtest;
            s.chgAh = obj.chgAh;
            s.disAh = obj.disAh;
            s.script = obj.script;
            s.temperature = obj.temperature;
            s.eta = obj.eta;
            s.QAh = obj.QAh;
            s.disZ = obj.disZ;
            s.disV = obj.disV;
            s.chgZ = obj.chgZ;
            s.chgV = obj.chgV;
        end
    end
end