classdef UnitTest < matlab.unittest.TestCase
% UnitTest is a modular test case class for MATLAB applications. 
% When executed via runtests(), this class will automatically run a series
% of unit tests against the specified executable, and create Markdown 
% reports summarizing the results. The MATLAB profiler is turned on during 
% unit test execution to quantify the code coverage.
%
% Author: Mark Geurts, mark.w.geurts@gmail.com
% Copyright (C) 2017 University of Wisconsin Board of Regents
%
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the  
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License along 
% with this program. If not, see http://www.gnu.org/licenses/.

    % Define UnitTest class properties
    properties

        % executable stores the application handle or directory name
        executable = 'fusion_tools'
        
        % configFile stores the name of the config file relative path
        configFile = ''
        
        % inputData stores a cell array of input datasets that will be
        % tested during these unit tests
        inputData = {
            'Head and neck'         '../test_data/CT HFS Pinnacle/CT/'
        }
    
        % referenceData stores a cell array of MATLAB variables containing
        % the expected contents of each archive above
        referenceData = {}
    
        % exportPath stores a temporary path to export files to
        exportPath = ''

        % figure stores the application UI handle. It is used during
        % teardown to close the figure
        figure = []
        
        % config stores a structure of config file contents
        config = struct
        
        % configOriginal stores a duplicate copy of the original config
        % file contents, to return the config file to its original state
        configOriginal = struct
        
        % reportFile stores the name of the markdown report file
        reportFile = '../test_results/unit_test_report.md'
        
        % traceMatrix stores the name of the markdown traceability matrix
        traceMatrix = '../test_results/unit_trace_matrix.md'
        
        % stats store the profiler results
        stats = []
    end
    
    % Define test level setup functions
    methods(TestMethodSetup)
        
    end
 
    % Define test level teardown functions
    methods(TestMethodTeardown)
        CloseFigure(testCase)
    end
 
    % Define class level setup functions
    methods (TestClassSetup)
        StartProfiler(testCase)
    end
    
    % Define class level teardown functions
    methods (TestClassTeardown)
        WriteMarkdownReport(testCase)
        WriteTraceabilityMatrix(testCase)
        StopProfiler(testCase)
    end
    
    % Define supporting methods
    methods (Static, Access = public)
        info = CPUInfo()
        info = MemInfo()
        sl = SLOC(file)
        image = LoadDICOMImages(path)
        varargout = StoreResults(varargin)
        modimage = ApplyRigidReg(image, v)
    end
    
    % Define unit tests
    methods(Test, TestTags = {'Unit'})
        Test01(testCase)
        Test02(testCase)
    end
end