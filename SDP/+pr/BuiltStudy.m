classdef BuiltStudy < com.TemperatureIndexedTestWrapper
    %BUILTSTUDY Stores data collected during a pulse-resistance study. 
    % 
    % A PR study could be conducted for a single cell at various 
    % temperatures, or several cells from the same production batch, each 
    % at a different temperature. 

    methods(Static)
        function merged = merge(name, varargin)
            %MERGE Combine two or more studies into a single study.
            %
            % MERGE(name,study1,study2,...) combines a set of PR studies
            % into one object with the specified name.

            testCntTotal = sum(cellfun(@(study)length(study.tests), varargin));
            tests = pr.Test.empty(0,testCntTotal);
            idxTest = 1;
            for idxStudy = 1:length(varargin)
                study = varargin{idxStudy};
                testCnt = length(study.tests);
                tests(idxTest:idxTest+testCnt-1) = study.tests;
                idxTest = idxTest + testCnt;
            end
            merged = pr.BuiltStudy(name, tests);
        end
    end
end