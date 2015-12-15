classdef regionBdd < baseRegion
%regionBdd is the punctured unit disk region.

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

properties
    dipole                      % Dipole source/sink location
end

methods
    function R = regionBdd(dv, qv, varargin)
        if ~nargin
            args = {};
        else
            args = [dv, qv, varargin];
        end
        R = R@baseRegion(args{:});
        if ~nargin
            return
        end
    end
end

methods(Access=protected)
    function m = mGetter(R)
        %Connectivity.
        m = numel(R.centers);
    end
    
    function sanityCheck(R)
        %Checks that region is in a valid state. Throws error if not.
        
        if ~isequal(numel(R.circulation), numel(R.centers))
            error(PoTk.ErrorTypeString.InvalidValue, ...
                ['Circulation vector must have the same number of\n'...
                'elements as there are inner boundaries.'])
        end
    end
end

end
