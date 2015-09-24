%% Potential class
% The <matlab:doc('./potential') |potential|> class handles all of the
% details of calculating the
% potential described in the <introduction.html introduction>.
% It operates in two stages.
% The first stage is invovled with precomputing the prime functions used for
% the potential, which depend solely on the definition of the
% <flowRegionDoc.html |flowRegion|>.
% The second stage is concerned with computing the values of the potential
% function at points in the flow domain. Using the potential function is
% then a two step process.

% Setup the required flow domain.
clear
islands = circleRegion({...
    circle(-2.03371+1.93258i, 0.92545), ...
    circle(3.78652+0.337079i, 0.876404), ...
    circle(-1.42697-2i, 0.908933)});
icirc = [-1, 0, 1];
vortices = [-3.2921-0.85393i, 0.75281+0.067416i, ...
    1.6517+3.2809i, 2.8427-2.6966i];
vcirc = [1, 1, 1, -1];
uniformFlow = 1;
uniformAngle = pi/4;
pd = flowRegion(islands, icirc, vortices, vcirc, ...
    uniformFlow, uniformAngle);

% Precomputes the required prime functions.
W = potential(pd);

% Now we can use it to evaluate points in the domain.
zeta = flowSamplePoints(pd, 25);
Wval = W(zeta);
velocity = velocityField(W, zeta);

disp([Wval(100); velocity(100)])


%%
% As shown in the <userguide.html quick summary>, the |potential| class also
% has builtin plot functionality. The default is to plot the stream lines,
% done by simply calling
%
%   plot(W)
%
% Equivalently one may use
%
%   plot(W, 'streamLines')
%
% or even
%
%   plotStreamLines(W)
%
% To plot a velocity vector field one may use
%
%   plot(W, 'velocityField')
%
% or 
% 
%   plotVelocityField(W)
