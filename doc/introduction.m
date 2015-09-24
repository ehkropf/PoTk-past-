%% PoTk Introduction
% This is a brief introduction to the theory that drives the
% Potential Toolkit (PoTk). For a much richer discussion, see
% 
% * list of appropriate references.

%% Complex potential flow
% We seek to compute a function
% \[ W(\zeta) = \varphi(\zeta) + i\psi(\zeta) \]
% where \(\varphi\) is a real-valued function known as the velocity potential
% and \(\psi\) is a real-valued function
% known as the stream function. By definition the flow velocity field is the
% gradient of this potential function, and in complex variables we can
% apply the Cauchy-Riemann equations, given \(\zeta = x + iy\), to write 
% \[ V(\zeta) = u(\zeta) + iv(\zeta) := \varphi_x + i\varphi_y = 
% \psi_y - i\psi_x. \]
%
% The potential computation can be split into two parts expressed by the
% linear combination
% \[ W(\zeta) = W_U(\zeta) + W_\Gamma(\zeta) \]
% representing, respectively, uniform background flow and flow due to
% \(n\) point vortices with circulation strengths
% \(\{\Gamma_k\}_{k=1}^n\).
% Let \(U\) be the strength of the background flow,
% let \(\chi\) be the angle of the
% background flow with respect to the positive real axis, and let
% \(\{\alpha_k\}_{k=1}^n\) be the locations of the point vortices.
%
% In the case there are no obstacles in the domain, these flows are given
% by
% \[ W_U(\zeta) = -iU\zeta e^{-i\chi} \]
% and
% \[ W_\Gamma(\zeta) = \sum_{k=1}^n \frac{\Gamma_k}{2\pi i} \log(\zeta -
% \alpha_k). \]
%
% If there are \(m>0\) distinct circular islands in the domain
% -- the domain is an \(m\)-connected, unbounded circle region --
% we must switch to the use of the
% Schottky-Klein prime function, \(\omega(\cdot,\cdot)\), to calculate the
% potential functions.

%% Unbounded, multiply connected domains
% Let \(\zeta(z)\) be a Mobius transformation which takes one of the
% \(m\) circles to the unit circle,
% the point at infinity to a point \(\beta\)
% inside the unit circle, and the rest of the boundary circles to circles
% inside the unit circle. Let \(\{\delta_j,q_j\}_{j=1}^{m-1}\) for
% be the centers and radii, respectively, of these
% reflected circles, and define the additional Mobius transformations
% \[ \theta_j(\zeta) := \delta_j + \frac{q_j^2\zeta}{1 -
% \overline{\delta_j}\zeta}. \]
% The points \(\alpha_k\) will now be images under
% \(\zeta\) of the \(n\) vortex locations in the \(z\)-plane,
% and we will allow for circulation around the circular boundaries, these
% circulation strengths given by \(\{\gamma_j\}_{j=0}^{m-1}\).
% Following [appropriate reference], we define the functions
% \[ G_0(\zeta,\alpha) := \frac{1}{2\pi i} \log\left(
% \frac{\omega(\zeta,\alpha)}{|\alpha|\omega(\zeta,1/\overline\alpha)}
% \right) \]
% and, for \(j\in \{1,\ldots,m-1\}\),
% \[ G_j(\zeta,\alpha) := \frac{1}{2\pi i} \log\left(
% \frac{\omega(\zeta,\alpha)}
% {|\alpha|\omega(\zeta,\theta_j(1/\overline\alpha))} \right). \]
%
% By [appropriate reference] we can then write down that
% \[ W_U(\zeta) = -4\pi U \left.\left( 
% \frac{\partial G_0(\zeta,\alpha)}{\partial \alpha_y} \cos\chi 
% + \frac{\partial G_0(\zeta,\alpha)}{\partial \alpha_x} \sin\chi
% \right)\right|_{\alpha=\beta}, \]
% where we have used \(\alpha = \alpha_x + i\alpha_y\), and
% \[ W_\Gamma(\zeta) = \sum_{k=1}^n \Gamma_k \left[ G_0(\zeta,\alpha_k)
% - G_0(\zeta,\beta) \right] - \sum_{j=0}^{m-1} \gamma_j G_j(\zeta,\beta). \]
%
% Numerically the \(G_j\) are never computed directly. A function theory
% identity is <nummethods.html instead employed> in combination with the
% already computed \(G_0\).

%% Bounded, multiply connected domains.
% These are different. Note
% \[ W_\Gamma(\zeta) = \sum_{k=1}^n \Gamma_k G_0(\zeta, \alpha_k) +
% \sum_{j=1}^m \gamma_j \left[ G_0(\zeta, \beta) - G_j(\zeta, \beta)
% \right]. \]
% This is due to putting all "extra" circulation on \(C_0\) instead of at
% \(\beta\).
