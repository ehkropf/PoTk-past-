%% Simple example
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
