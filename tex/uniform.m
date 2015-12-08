%% Uniform flow examples.


%%
% Some setup for use below.

clear
aspectequal = @(ax) set(ax, 'dataaspectratio', [1, 1, 1]);
streamColor = [0, 0.447, 0.741];
vectorColor = [0.929, 0.694, 0.125];


%% The most basic case.
% The vector components of uniform flow at speed \(U\)in the plane at an
% angle \(\chi\) with respect to the real axis are given by
% \[ u = U\cos\chi \quad\mathrm{and}\quad v = U\sin\chi, \]
% from which we get
% \[ \frac{dW_U}{dz} = Ue^{-i\chi} \]
% so that the complex potential for the uniform flow is
% \[ W_U(z) = Uze^{-i\chi}. \]
% Let's view a flow 

U = 1;
chi = pi/4;

xylim = [-1, 1, -1, 1];
z = rectgridz(xylim, 200);
zv = rectgridz(xylim, 20, 0.01);

Wz = U*exp(-1i*chi)*z;
vzv = 1i*U*exp(-1i*chi)*ones(size(zv));

quiver(real(zv), imag(zv), real(vzv), imag(vzv), 'color', vectorColor)
hold on
contour(real(z), imag(z), imag(Wz), 20, 'lineColor', streamColor)
hold off
aspectequal(gca)
axis(xylim)
axis off


%% Flow past one object.
% Uniform flow past a circle is given by the Milne-Thompson circle theorem
% via
% \[ W(z) = U\left( z + \frac{1}{z} \right). \]

U = 1;
Wu = @(z, U) U*(z + 1./z);
dWdz = centdiffz(@(z) Wu(z, U));

xylim = 2.5*xylim;
z = rectgridz(xylim, 200);
outer = abs(z) > 1;
zv = rectgridz(xylim, 20, 0.01);
outv = abs(zv) > 1;

Wz = complex(nan(size(z)));
Wz(outer) = Wu(z(outer), U);
vzv = complex(nan(size(zv)));
vzv(outv) = dWdz(zv(outv));

quiver(real(zv), imag(zv), real(vzv), imag(vzv), 'color', vectorColor)
hold on
contour(real(z), imag(z), imag(Wz), 20, 'lineColor', streamColor)
fill(circle(0, 1))
plot(circle(0, 1))
hold off
aspectequal(gca)
axis(xylim)


%%
% This may also be modified by a line vortex via
% \[ W(z) = U\left( z + \frac{1}{z} \right) + \frac{\Gamma}{2\pi i}
% \log z \]
% since the effect is to maintain the boundary of the unit circle as a
% streamline.

Gamma = -3;
Wc = @(z, Gamma) Gamma/(2i*pi)*log(z);
W = @(z, U, Gamma) Wu(z, U) + Wc(z, Gamma);
dWdz = centdiffz(@(z) W(z, U, Gamma));

Wz(outer) = W(z(outer), U, Gamma);
vzv(outv) = dWdz(zv(outv));

quiver(real(zv), imag(zv), real(vzv), imag(vzv), 'color', vectorColor)
hold on
contour(real(z), imag(z), imag(Wz), 20, 'lineColor', streamColor)
fill(circle(0, 1))
plot(circle(0, 1))
hold off
aspectequal(gca)
axis(xylim)


%%
% Stagnation point at \(z = e^{i3\pi/2}\).

Gamma = -4*pi;
dWdz = centdiffz(@(z) W(z, U, Gamma));

Wz(outer) = W(z(outer), U, Gamma);
vzv(outv) = dWdz(zv(outv));

quiver(real(zv), imag(zv), real(vzv), imag(vzv), 'color', vectorColor)
hold on
contour(real(z), imag(z), imag(Wz), 20, 'lineColor', streamColor)
fill(circle(0, 1))
plot(circle(0, 1))
hold off
aspectequal(gca)
axis(xylim)


%%
% Stagnation point away from the unit disk.

Gamma = -5*pi;
dWdz = centdiffz(@(z) W(z, U, Gamma));

xylim = [-2.5, 2.5, -3.5, 1.5];
z = rectgridz(xylim, 200);
outer = abs(z) > 1;
zv = rectgridz(xylim, 20, 0.01);
outvec = abs(zv) > 1;

Wz = complex(nan(size(z)));
Wz(outer) = W(z(outer), U, Gamma);
vzv = complex(nan(size(zv)));
vzv(outvec) = dWdz(zv(outvec));

quiver(real(zv), imag(zv), real(vzv), imag(vzv), 'color', vectorColor)
hold on
contour(real(z), imag(z), imag(Wz), 20, 'lineColor', streamColor)
fill(circle(0, 1))
plot(circle(0, 1))
hold off
aspectequal(gca)
axis([-2.5, 2.5, -3.5, 1.5])
