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
        
        sanityCheck(R)
    end
end

methods(Access=protected)
    function m = mGetter(R)
        m = numel(R.centers);
    end
    
    function sanityCheck(R)
        % Add checks to base set.
        
        sanityCheck@baseRegion(R)
    end
end

end
