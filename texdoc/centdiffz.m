function dw = centdiffz(f, h)
%centdiffz is the central difference formula for f(z).

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

if nargin < 2
    h = 1e-8;
end

dh = 0.5*h;
dw = @(z) (imag(f(z + 1i*dh) - f(z - 1i*dh)) ...
    - 1i*imag(f(z + dh) - f(z - dh)))/h;
