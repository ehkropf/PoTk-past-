%% PoTk Numerical Methods.
% Various methods for computational shortcuts are described here.


%% Definitions.
% We recall that
% \[ G_0(\zeta,\alpha) := \frac{1}{2\pi i} \log\left(
% \frac{\omega(\zeta,\alpha)}{|\alpha|\omega(\zeta,1/\overline\alpha)}
% \right) \]
% and
% \[ G_j(\zeta,\alpha) := \frac{1}{2\pi i} \log\left(
% \frac{\omega(\zeta,\alpha)}
% {|\alpha|\omega(\zeta,\theta_j(1/\overline\alpha))} \right). \]
% Also recall that there is a unique function \(X(\zeta,\alpha)\) such that
% \[ \omega(\zeta,\alpha) = \sqrt{X(\zeta,\alpha)}. \]
% The branch of the square root is chosen so that \(\omega(\zeta,\alpha)\)
% behaves like \((\zeta - \alpha)\) as \(\zeta\to\alpha\).
%
% A defining property of \(X\) is the Hejhal identity [reference]
% \[ X(\theta_j(\alpha), \zeta) = \exp\left( -4\pi i \left[ v_j(\alpha)
% - v_j(\zeta) + \tfrac{1}{2} \tau_{jj} \right] \right) \theta_j'(\alpha)
% X(\alpha,\zeta). \]
% We recall also that \(X(\zeta,\alpha) = X(\alpha,\zeta)\), so in taking the
% square root we find
% \[ \omega(\zeta, \theta_j(\alpha)) = -\omega(\zeta, \alpha) \exp\left(
% -2\pi i \left[ v_j(\alpha) - v_j(\zeta) + \tfrac{1}{2} \tau_{jj} \right]
% \right) \sqrt{\theta_j'(\alpha)} \]
% where the negative sign was chosen to ensure the proper branch of the
% square root.
%
% Given a point \(z\in\{\mathbb{C}\cup\infty\}\), define its conjugate
% point with respect to the unit cirlce to be
% \[ z^* := 1/\overline{z}. \]


%% Computing \(G_j\).
% Applying the Hejhal identity we can deduce that
% \[ G_j(\zeta, \alpha) = \frac{1}{2\pi i} \log\left( \frac{\omega(\zeta,
% \alpha)}{|\alpha|\omega(\zeta, \alpha^*)} \cdot 
% \frac{-exp\left( 2\pi i \left[ v_j(\alpha^*) - v_j(\zeta) + \tfrac{1}{2}
% \tau_{jj} \right] \right)}{\sqrt{\theta_j'(\alpha^*)}} \right). \]
% This form is used instead of the more obvious \(G_0(z,\alpha) + v_j(z) +
% \cdots\), since numerically there is no branch cut confusion with the
% given formula.
