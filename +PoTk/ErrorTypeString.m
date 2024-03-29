classdef ErrorTypeString
%PoTk.ErrorTypeString encapsulates error ID strings.

% Everett Kropf, 2015
% 
% This file is part of the Potential Toolkit (PoTk).
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

% properties(Dependent)
%     InvalidArgument
% end

methods(Static)
    function str = InvalidArgument()
        str = 'PoTk:InvalidArgument';
    end
    
    function str = InvalidValue()
        str = 'PoTk:InvalidValue';
    end
    
    function str = RuntimeError()
        str = 'PoTk:RuntimeError';
    end
    
    function str = UndefinedState()
        str = 'PoTk:UndefinedState';
    end
end

end
