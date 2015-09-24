%% PoTk User's Guide
% Welcome to the user's guide for the Potential toolkit (PoTk) for MATLAB.
%
% * <#1 Overview>
% * <introduction.html Introduction>
% * <componentList.html Component list>
%
% This page presents a quick summary regarding the use of the toolkit.
% See the the documentation on the individual commands
% for more information than may be found
% here.

%% Overview
% Use of the toolkit generally falls into the following pattern:
%
% # Defining the flow domain via use of the |flowRegion| class.
% # Precomputing the potential function by calling the |potential| class
% constructor with the flow domain.
% # Either evaluating the potential function at given points or getting
% the flow field values at given points.
%
% -- OR --
%
% # Using |potool| (the PoTk Interactive Figure), which handles the
% above steps automatically.

%% Flow domain
% The |flowRegion| class describes the physical domain for the problem.
% A |circleRegion| object is used to describe any circular islands
% (circular obstacles; regions of no flow) in the domain,
% and is built using |circle| objects. The |circleRegion| and |circle|
% classes are provided via the Conformal Mapping Toolbox (CMT).

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
plot(pd)


%% Potential
% Once the flow domain is defined, we can setup the potential function.
% This involves calculating the Schottky-Klien prime functions
% \(\omega(\cdot,\cdot)\) used in the representation of the potential.

W = potential(pd);


%%
% *Evaluation*
%
% As an example of evaluating the potential, at points in the plane,
% we can ask for a plot of the stream lines.

plot(W)


%%
% We can also ask for a plot of the
% velocity vector field. Note the stream lines are
% also plotted. This can be controlled by various options.

plot(W, 'velocityField')


%%
% Individual points (or vectors of points) may also be evaluated, both for
% the value of the potential function, or the velocity vector values.

disp(W(0))

%%
disp(velocityField(W, 0))


%% Graphical user interface
% All of the above steps, with the exception of evaluating the potential at
% arbitrary points, may be automated by using the interactive figure.
% Simply type
%
%   potool
%
% at the command line to launch the tool.
%
% <<../potool_image.png>>
%
% The user may
% 
% * Click in the axes to start adding a circular obstacle. Moving the mouse
% changes the circle radius, and a second click finishes adding the circle.
% Pressing the |escape| key during this process cancels adding the circle.
% * Shift+click in the axes to add a point vortex.
% * Press the *Axis* buton to change the axes limits.
% * Click the *Reset* button to clear the axes.
% * Click the *Delete* button to delete a circle or vortex. Clicking near a
% circle or vortex completes this process, pressing the |escape| key
% cancels the delete.
% * Click the *Export* button to export variables from the GUI to the current
% workspace.
% * Click the *Stream lines* button to calculate the potential and plot
% the stream lines.
% * Click the *Velocity field* button to calculate the potential, take
% its derivative, and plot the velocity vector field. 
% This may be done with or without the stream lines.
% * Use the _Island_ and _Vortex_ panels to change the scalar which
% dictates the circulation around obstacle circles and point vortices
% respectively. The circle location and radii may be edited, as may the
% vortex locations.
% * Use _Background flow_ panel to set the strength and angle of the uniform
% background flow.
