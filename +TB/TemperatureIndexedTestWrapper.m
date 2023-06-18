classdef TemperatureIndexedTestWrapper < handle
    %TEMPERATUREINDEXEDTESTWRAPPER Generic wrapper for tests (and similar) 
    % objects that correspond to temperatures. Often used as a base class
    % for Study-type objects.

    properties
        name             % Name of the set.
        testType         % Test class name.
        tests            % Vector of test (or similar) objects.
        testTemperatures % Temperatures corresponding to each test [C].
        testNames        % Name of each test (for simpler lookup).
    end

    methods
        function obj = TemperatureIndexedTestWrapper(name, tests)
            obj.name = name;
            obj.tests = tests;
            obj.testType = class(tests);

            % Sort tests by temperature.
            [~, idx] = sort([tests.temp]);
            tests = tests(idx);
            obj.tests = tests;

            obj.testTemperatures = [tests.temp];
            obj.testNames = {tests.name};

            % Check for duplicate temperatures.
            if length(obj.testTemperatures) ~= length(unique(obj.testTemperatures))
                warning(['Duplicate tests (at same temperature) ' ...
                    'found for %s.'], name);
            end
        end

        function tests = getTest(obj, T, name)
            %GETTEST Fetch the test(s) performed at the specified
            %  temperature and name or raise an error if there is 
            %  no matching test.
            %
            % tests = study.GETTESTS(T) fetches all tests at temperature T.
            %
            % tests = study.GETTESTS(T,name) fetches all tests at
            %   temperature T and named NAME.

            if nargin < 3
                tests = obj.tests(obj.testTemperatures==T);
            else
                name = string(name);
                tests = obj.tests(and(obj.testTemperatures==T,obj.testNames==name));
            end

            if isempty(tests)
                error('Tests(s) with matching temperature/name not found.');
            end
        end % getTest()

        function s = toStruct(obj)
            s.name = obj.name;
            s.testTemperatures = obj.testTemperatures;
            s.testNames = obj.testNames;
            for k = length(obj.tests):-1:1
                t(k) = obj.tests(k).toStruct();
            end
            s.tests = t;
        end
    end
end