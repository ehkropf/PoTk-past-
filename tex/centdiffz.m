function dw = centdiffz(f, h)
%centdiffz is the central difference formula for f(z).

% Copyright Everett Kropf, 2015.

if nargin < 2
    h = 1e-8;
end

dh = 0.5*h;
dw = @(z) (imag(f(z + 1i*dh) - f(z - 1i*dh)) ...
    - 1i*imag(f(z + dh) - f(z - dh)))/h;
