%% The modified Green's function, the Schottky-Klein prime function, and Mobius transformations
% We'll start with an unbounded circle domain, which will be denoted \(D\).
% Note that in the code used,
% circles are represented by a |circle| object which is constructed by
% the syntax 
%
%   acirc = circle(center, radius);
%
% (The circles, along with the circle region are part of the Conformal Mapping
% Toolbox (CMT) which may be found
% <https://github.com/tobydriscoll/conformalmapping here>.)
% We'll put in a point vortex given by the |vortex| variable, and make some
% test points. A picture is given with the vortex represented by a blue 'X', and
% the test points as black dots. The origin is shown as a red '+'.

clear
islands = circleRegion(...
    circle(-2.1757+1.23364i, 1.25387), ...
    circle(1.09907-1.97383i, 1.47459), ...
    circle(3.20748+3.2972i, 1.34748));
s = islands.centers;
r = islands.radii;
m = numel(s);

vortex = 0.958884414225943 + 1.11622384937239i;
testPts = [-2.06840905857741 - 1.91106962343096i
          -0.727460271966526 + 4.5295480334728i
            4.51443043933054 - 1.15932560669456i
         -0.0366684728033473 + 0.283210209205022i
            1.89348508368201 + 1.6241589958159i];

invfill(islands)
hold on
for j = 1:m
    text(real(s(j)), imag(s(j)), sprintf('%d', j-1))
end
plot(islands)
plot(0, 0, 'r+')
plot(vortex, 'bx')
plot(testPts, 'k.')
hold off

       
%%
% We'll now use an affine map to make the first circle unit, and then
% compose this with a Mobius transformation which puts the other circles
% inside the first. The composition map will be called \(\zeta_1(z)\). The
% image of the island domain will be labeled \(D_1\).

zunit = mobius(1, -s(1), 0, r(1));
zeta1 = mobius(0, 1, 1, 0)*zunit;
display(zeta1)

D1 = zeta1(islands);
alph1 = zeta1(vortex);
z1test = zeta1(testPts);

invfill(D1)
hold on
for j = 2:m
    text(real(D1.centers(j)), imag(D1.centers(j)), sprintf('%d', j-1))
end
plot(D1)
plot(alph1, 'bx')
plot(z1test, 'k.')
hold off


%%
% We'll now define the modified Green's function with respect to the first
% circle in the usual way via the prime function and a Mobius transformation
% which takes domain \(D\) to a region bounded by the unit circle.
% That is,
% \[
%  G_0(\zeta, \gamma) = \frac{1}{2\pi i} \log\left(
%  \frac{\omega(\zeta, \gamma)}{|\gamma|
%  \omega(\zeta, 1/\overline{\gamma})} \right).
% \]
% This is done using an anonymous function, where, due to the way the prime
% function is calculated below, we'll use the functions |wt| to represent
% the prime function in the numerator and |wb| for the prime function in
% the denominator.

G0 = @(z, gamma, wt, wb) log(wt(z)./wb(z)/abs(gamma))/(2i*pi);


%%
% Now we use the \(\zeta_1(z)\) map to calculate the value of this Green's
% function at the test points.

omeg1 = skprime(alph1, D1);
omeg1o = invParam(omeg1);
wz1 = G0(z1test, alph1, omeg1, omeg1o);


%%
% We now compute a Mobius transformation, \(m_{ann}\) which takes the second
% disk in \(D_1\) to the inner boundary of an annulus by finding a value
% \(a\) such that
% \[
%  m_{ann}(z) := \frac{z - a}{1 - z\overline{a}}.
% \]
% The resulting composition
% \[
%  \zeta_2(z) := \left( m_{ann} \circ \zeta_1 \right) (z)
% \]
% is then a map from the unbounded region \(D\) to an annular region which
% we will denote \(D_2\).

d2 = D1.centers(2);
q2 = D1.radii(2)^2;
k = exp(1i*angle(d2));
d2 = abs(d2)^2;
a = k*(1 + d2 - q2 - sqrt((d2 - q2)^2 - 2*(d2 + q2) + 1))/(2*sqrt(d2));

mann = mobius(1, -a, -conj(a), 1);
display(mann)
zeta2 = mann*zeta1;
display(zeta2)

D2 = zeta2(islands);
alph2 = zeta2(vortex);
z2test = zeta2(testPts);

invfill(D2)
hold on
for j = 2:m
    text(real(D2.centers(j)), imag(D2.centers(j)), sprintf('%d', j-1))
end
plot(D2)
plot(alph2, 'bx')
plot(z2test, 'k.')
hold off


%%
% Now we use the \(\zeta_2\) map to compute the associated modified Green's
% function.

omeg2 = skprime(alph2, D2);
omeg2o = invParam(omeg2);
wz2 = G0(z2test, alph2, omeg2, omeg2o);


%%
% Comparing the results shows that, with the two domains which differ by a
% Mobius transformation, the values of the modified Green's function differ
% only by a real constant. We'd like the real constant to be zero, but this
% is the next best thing.

display(num2str(wz2 - wz1, '%.5f'))


%%

max(diff((real(wz2 - wz1))))
norm(imag(wz2 - wz1), inf)
