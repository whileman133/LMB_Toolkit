classdef BuiltStudy < TB.TemperatureIndexedTestWrapper
    %BUILTSTUDY Stores data collected during an OCP study. 
    % 
    % An OCP study could be conducted for a single cell at various 
    % temperatures, or several cells from the same production batch, each 
    % at a different temperature. 
    %
    % Each cell dataset is called a 'test'. Provide a test at the reference
    % temperature (temperatureRef) to allow calibrating coloumbic
    % efficiency and calculating relative capacity.

    properties(SetAccess=protected)
        temperatureRef   % Calibration reference temperature, usu. 25 [C].
        crate            % C-rate at which the study was performed.
    end

    methods
        function obj = BuiltStudy(name, Tref, crate, tests)
            %OCPSTUDY Construct an OCP study from component tests.

            obj = obj@TB.TemperatureIndexedTestWrapper(name,tests);
            obj.temperatureRef = Tref;
            obj.crate = crate;

            temperatures = obj.testTemperatures;
            tests = obj.tests;

            % Check for missing test at calibration temprature.
            if ~any(temperatures == Tref)
                warning(['Missing test at reference temperature T=%d. ' ...
                    'Skipping calibration for %s'], Tref, name);
                return;
            end

            % Calibration -------------------------------------------------

            % First determine eta and Q for test(s) at T=Tref.
            testsTref = tests(temperatures == Tref);
            etaTref = zeros(size(testsTref));
            for k = 1:length(testsTref)
                testTref = testsTref(k);
                s12DisAh = testTref.disAh(testTref.script<=2);
                s12ChgAh = testTref.chgAh(testTref.script<=2);
                etaTref(k) = testTref.disAh(end) / testTref.chgAh(end);
                Q = s12DisAh(end) - etaTref(k) * s12ChgAh(end);
                testTref.processOCP(etaTref(k), Q);
            end

            % Use average eta(Tref) to process remaining tests.
            etaTref = mean(etaTref);

            % Now determine eta and Q for remaining tests.
            for test = tests(temperatures ~= Tref)
                s12DisAh = test.disAh(test.script<=2); % Ah discharged up to scripts 1 and 2.
                s12ChgAhTref = test.chgAhTref(test.script<=2); % Ah charged up to scripts 1 and 2 at Tref.
                s12ChgAhTtest = test.chgAhTtest(test.script<=2); % Ah charged up to scripts 1 and 2 at Ttest.
                eta = (test.disAh(end) - etaTref * test.chgAhTref(end)) / test.chgAhTtest(end);
                Q = s12DisAh(end) - etaTref * s12ChgAhTref(end) - eta * s12ChgAhTtest(end);
                test.processOCP(eta, Q);
            end
        end
    end
end