%% PoTk: a Guide of Sorts
% Welcome to some sort of guide for the Potential Toolkit (PoTk) for MATLAB.


%% Introduction
% There should be some introductory material here. Maybe we explain the need
% for the toolkit, or say some things about potential flow or potential
% problems in the complex plane.
%
% For a test of documentation setup, we need equations like \(y = mx + b\),
% or some silly thing. We'll also test with a separated equation.
% \[ \log(x) := \int_1^x \frac{1}{t} \,dt. \]
%
% Note the PoTk relies heavily on the Conformal Mapping Toolkit
% (<http://github.com/tobydriscoll/conformalmapping CMT>) by Everett Kropf
% and Toby Driscoll. It will be mentioned often in the subsequent.


%% The physical domain
% All potential calculations are done over a physical flow domain. The
% MATLAB class |flowRegion| encapsulates this. By the nature of the
% toolkit, any obstacles in the flow domain are circles. Thus the basic
% unit of a boundary is the |circle| object from the CMT. We will be
% interested in multiple circles in general, thus the first object to
% consider is a |circleRegion| class, also defined in the CMT.

clear
islands = circleRegion({...
    circle(-1.98876+1.79775i, 1.696), ...
    circle(2.68539-1.30337i, 1.66505), ...
    circle(-1.98876-2.74157i, 0.902521), ...
    circle(1.85393+1.88764i, 0.660919)});
plot(islands)


%%
% A basic flow domain exterior to the above 4 circles as obstacles would
% then be given by the following. Note we have specified a circulation
% strength for every island defined. This is required by |flowRegion|.

pd = flowRegion(islands, [-1, 1, 0, -2]);
plot(pd)


%%
% Point vortices and their circulation strengths may also be specified.

clear
islands = circleRegion({...
    circle(-2.32584+1.95506i, 1.42834), ...
    circle(3.60674+0.314607i, 1.25782), ...
    circle(-2.16854-2.5618i, 1.02586)});
icirc = [0, 0, 0];
vortices = [1.2247-3.0112i, -4.0337-1.0787i, 1.2697+2.1348i];
vcirc = [1, -1, 1];
pd = flowRegion(islands, icirc, vortices, vcirc);
plot(pd)


%%
% Uniform flow strength and angle with respect to the \(x\)-axis may also
% be given.

U = 1;
Chi = pi/4;
pd2 = flowRegion(islands, icirc, vortices, vcirc, U, Chi);


%%
% As given in the examples above, these arguments are positional for the
% |flowRegion| constructor. They can, however, also be given in the
% name/value pair argument form.

pd3 = flowRegion(...
    'islands', islands, ...
    'islandCirculation', icirc, ...
    'vortexLocation', vortices, ...
    'vortexCirculation', vcirc, ...
    'uniformStrength', U, ...
    'uniformAngle', Chi);


%% Point vortex example
% Given a set of \(n\) vortices located at \(\alpha_k\), for
% \(k=\{1,\ldots,n\}\), with strength \(\Gamma_k\),
% the potential function for no obstacles in the flow is given by
% \[ W(\zeta) = \sum_{k=1}^n \frac{\Gamma_k}{2\pi i} \log \left( \zeta -
% \alpha_k \right). \]
% Given a flow domain with three point vortices, we calculate the potential
% and plot the streamlines.

clear
vortices = [1.2247-3.0112i, -4.0337-1.0787i, 1.2697+2.1348i];
Gamma = [1, -1, 1];
pd = flowRegion(...
    'vortexLocation', vortices, ...
    'vortexCirculation', Gamma);
W = potential(pd);

clf
plot(W)


%%
% It is also possible to plot the velocity vector field.

clf
plot(W, 'velocityField')


%% Point vortex outside a cylinder
% Now we compute the velocity field for a point vortex with unit
% circulation outside a cylinder (circular island).
% Without loss of generality (via a Mobius tansformation), we can
% suppose that the cylinder is the unit circle, and we let
% \[ \zeta(z) := \frac{1}{z}, \]
% and let \(\alpha\) be the image of the point vortex location under
% \(\zeta\).
% It is shown in [appropriate reference] that the potential is then
% given by
% \[ W(\zeta(z)) = \frac{1}{2\pi i} \log\left( 
% \frac{\zeta - \alpha}{|\alpha|(\zeta - 1/\overline\alpha)} \right)
% + \frac{1}{2\pi i} \log\zeta + K \]
% where \(K\) is some constant which we can take to be zero.

clear
pd = flowRegion(...
    circleRegion(circle(1.69663-0.404494i, 1.1443)), 0, ...
    -1.8764+1.8876i, 1);
W = potential(pd);

clf
plot(W, 'velocityField')


%% Multiple point vortices and cylinders
% By the use of the Shottky-Klein prime function, \(\omega(\cdot,\cdot)\), 
% we can write down the
% complex potential for multiple cylinders, with or without circulation,
% and multiple point vortices.
%
% Let now \(\zeta(z)\) be a Mobius transformation which takes one of the
% \(m\) cylinders to the unit circle,
% the point at infinity to a point \(\beta\)
% inside the unit circle, and the rest of the boundary circles to circles
% inside the unit circle. Let \(\delta_j,q_j\) for \(j\in
% \{1,\ldots,m-1\}\) be the centers and radii, respectively, of these circles,
% and define the Mobius transformations
% \[ \theta_j(\zeta) := \delta_j + \frac{q_j^2\zeta}{1 -
% \overline{\delta_j}\zeta}. \]
% The points \(\alpha_k\) will be images under
% \(\zeta\) of the \(n\) vortex locations in the \(z\)-plane.
%
% Following [appropriate reference], we define
% \[ G_0(\zeta,\cdot) := \frac{1}{2\pi i} \log\left(
% \frac{\omega(\zeta,\cdot)}{|\cdot|\omega(\zeta,1/\overline\cdot)}
% \right) \]
% and
% \[ G_j(\zeta,\cdot) := \frac{1}{2\pi i} \log\left(
% \frac{\omega(\zeta,\cdot)}{|\cdot|\omega(\zeta,\theta_j(1/\overline\cdot))}
% \right). \]
% Then for vortex circulation strengths \(\Gamma_k\) and cylinder
% circulation strengths \(\gamma_j\), we have the potential function
% \[ W(\zeta) = \sum_{k=1}^n \Gamma_k \left( G_0(\zeta,\alpha_k)
% - G_0(\zeta,\beta) \right) - \sum_{j=1}^{m-1} \gamma_j G_j(\zeta,\beta). \]

clear
islands = circleRegion({...
    circle(-2.03371+1.93258i, 0.92545), ...
    circle(3.78652+0.337079i, 0.876404), ...
    circle(-1.42697-2i, 0.908933)});
icirc = [-1, 0, 1];
vortices = [-3.2921-0.85393i, 0.75281+0.067416i, ...
    1.6517+3.2809i, 2.8427-2.6966i];
vcirc = [1, 1, 1, -1];

W = potential(flowRegion(islands, icirc, vortices, vcirc));
plot(W, 'velocityField')


%% Uniform flow
% According to [appropriate reference] the complex potential in the
% presence of a background uniform flow may be written
% \[ W = W_U + W_\Gamma \]
% where \(W_U\) is the potential for the uniform flow and \(W_\Gamma\) is
% the potential due to point vortex and cylinder circulation.
%
% In the case of no obstacles in the flow we may simply write
% \[ W_U(\zeta) = -iU\zeta e^{-i\chi} \]
% where \(U\) is the flow strength and \(\chi\) is the angle of the flow
% with respect to the posistive real axis.
%
% In the case of \(m>0\) cylinders, we resort to the prime function to
% calculate uniform flow. For \(G_0(\zeta,\alpha)\) we denote \(\alpha =
% \alpha_x + i\alpha_y\). Given that \(\zeta(z)\) is the Mobius
% transformation described above with \(\beta = \zeta(\infty)\) it can be
% shown that
% \[ W_U(\zeta) = -4\pi U \left.\left( \frac{\partial G_0}{\partial \alpha_y}
% \cos\chi - \frac{\partial G_0}{\partial \alpha_x} \sin\chi
% \right)\right|_{\alpha=\beta} \]
% where \(G_0(\cdot,\cdot)\) is as defined above.

U = 1;
chi = pi/4;
pd = flowRegion(islands, icirc, vortices, vcirc, U, chi);
W = potential(pd);
plot(W, 'velocityField')
