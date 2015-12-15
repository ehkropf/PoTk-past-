function tf = isRealAndFinite(v)
%isRealAndFinite verifies numeric, finite, and not-imaginary properties.
%
% tf = isRealAndFinite(v) is equivalent to calling
%
%    tf = PoTk.isFiniteNumber(v) & imag(v) == 0;
%
% See also PoTk.isFiniteNumber, imag.

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

tf = PoTk.isFiniteNumber(v) & imag(v) == 0;
