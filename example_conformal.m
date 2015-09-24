%% cmt/potk technology preview
clear

streamColor = [0, 0.447, 0.741];
vectorColor = [0.929, 0.694, 0.125];


%%
% Going to map the flow domain with vortices to flow around some shape via
% conformal mapping. First we need a domain, a potential, and a conformal map.

Df = flowRegion(...
    circleRegion(circle(0, 1)), 2, ...
    [1.5+0.5i, 1.5-0.5i], 3*[-1, 1], ...
    1, 0);

W = potential(Df, 'useWaitBar', true);

G = splinep([ ...
    1.5497 + 0.0058i; 1.5146 + 0.5088i; 0.7661 + 0.4971i
    -1.4795 - 0.0058i; 0.7661 - 0.4854i; 1.5263 - 0.4854i ...
]);
g = szmap(G)';


%%
% Get the sample points and values to plot this stuff.

xylim = plotbox(Df, 1.5);
np = 400;
% nv = 20;

zs = flowSamplePoints(Df, np, xylim, 0.01);
sl = complex(nan(size(zs)));
sl(~isnan(zs)) = W(zs(~isnan(zs)));

% zv = flowSamplePoints(Df, nv, xylim, 0.05);
% V = velocityField(W, zv);


%%
% Plot the circle domain to check.

figure(10), clf

subplot(1, 2, 1)
hold on
% quiver(real(zv), imag(zv), real(V), imag(V), 'color', vectorColor)
contour(real(zs), imag(zs), imag(sl), 20, 'linecolor', streamColor)
plot(Df)
axis(xylim)
aspectequal
axis off


%%
% Transplant to the conformal image domain and plot.

ws = g(zs);
wvort = g(Df.vortexLocation);

subplot(1, 2, 2)
hold on
axis(cmt.plotbox(ws, 1))
contour(real(ws), imag(ws), imag(sl), 20, 'linecolor', streamColor)
fill(G)
plot(G)
plot(real(wvort), imag(wvort), 'k.', 'markersize', 18)
aspectequal
axis off
