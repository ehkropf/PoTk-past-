function pobj = potential(region, varargin)
%POTENTIAL generates a complex potential object.
%
% P = potential(region, ...)
% Create complex potential object based on region type. If region is not
% of class |regionExt| or |regionBdd|, then the builtin MATLAB potential
% function is called.
%
% See also regionExt, potentialExt, regionBdd, potentialBdd, potential.

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

switch class(region)
    case 'regionBdd'
        pobj = potentialBdd(region, varargin{:});
        
    case 'regionExt'
        pobj = potentialExt(region, varargin{:});
        
    otherwise
        error(PoTk.ErrorTypeString.RuntimeError, ...
            'Unrecognized region type ''%s'' for potential.', ...
            class(region))
end
