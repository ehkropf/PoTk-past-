function a = annulusConst(d, q)
%annulusConst gives the self-map to annulus constant.

% E. Kropf, 2015
% 
% This file is part of SKPrime.
% 
% SKPrime is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% SKPrime is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with SKPrime.  If not, see <http://www.gnu.org/licenses/>.

% Save circle angle.
k = exp(1i*angle(d));

% Effective rotation of C_j to positive real axis.
d2 = abs(d)^2;
q2 = q^2;

a = (1 + d2 - q2 - sqrt((d2 - q2)^2 - 2*(d2 + q2) + 1))/(2*sqrt(d2));
assert(a > 0, 'Unexpected non-positive value detected.')
a = k*a;
