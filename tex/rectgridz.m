function z = rectgridz(xylim, res, pad)
%rectgridz contructs rectangular complex grid.

% E. Kropf, 2015

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
