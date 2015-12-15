%% Simple example

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

clear

D = circleRegion(...
    circle(-1.5028+1.83925i, 1.65815), ...
    circle(1.68224-2.31028i, 1.3772));
icirc = [-1, 1];
vorts = [3.1850+2.4224i; -5.0019+2.8935i; -3.6785-2.7364i];
vcirc = [1, -1, 1];
U = 0.08;

Df = flowRegion(D, icirc, vorts, vcirc, U);

W = potential(Df);

clf
plot(W, 'velocity')
