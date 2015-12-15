function z = rectgridz(xylim, res, pad)
%rectgridz contructs rectangular complex grid.

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
    res = 200;
end
if numel(res) == 1
    resx = res;
    resy = res;
else
    resx = res(1);
    resy = res(2);
end

if nargin < 3
    pad = 0;
end
if pad ~= 0
    padx = pad*diff(xylim(1:2));
    pady = pad*diff(xylim(3:4));
    xylim = [xylim(1) + padx, xylim(2) - padx, ...
        xylim(3) + pady, xylim(4) - pady];
end

[x, y] = meshgrid(linspace(xylim(1), xylim(2), resx), ...
    linspace(xylim(3), xylim(4), resy));
z = complex(x, y);
