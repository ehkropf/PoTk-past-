function result = potests(select)
%potests runs unit tests for PoTk.
%
% potests
%   Runs all the available tests for the prime function in the test suite.
%
% potests <string>
%   Runs a test in the test suite specified by <string>. Do not prepend the
%   string with 'PoTkUnitTest.' as this is done automatically. The wildcard
%   character is accepted. See the TestSuite.fromPackage documentation
%   for more information on string format.
%
% See also: matlab.unittest.TestSuite.fromPackage

% Copyright Everett Kropf, 2015
% 
% This file is part of PoTk.
% 
% PoTk is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% PoTk is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with PoTk.  If not, see <http://www.gnu.org/licenses/>.

import matlab.unittest.TestSuite
import matlab.unittest.TestRunner

runner = TestRunner.withTextOutput('Verbosity', 2);

suiteArgs = {'Name', 'PoTkUnitTest.*'};
if nargin
%     switch select
%         case 'prime'
%             % Now set as default.
%         case {'g0', 'G0'}
%             suiteArgs = {'Name', 'PoTkUnitTest.G0*'};
%         case 'all'
%             suiteArgs = {};
%         otherwise
            suiteArgs = {'Name', ['PoTkUnitTest.' select]};
%     end
end
tests = TestSuite.fromPackage('PoTkUnitTest', suiteArgs{:});

rng('shuffle')
result = run(runner, tests);

fprintf('\nTest run summary:\n\n')
disp(table(result))

if ~nargout
    clear result
end
