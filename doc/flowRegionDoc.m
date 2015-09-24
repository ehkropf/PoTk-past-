%% Flow domain class
% The <matlab:doc('flowRegion') |flowRegion|> is defined by zero or more of
% the following:
%
% * An unbounded circle region composed of zero or more non-intersecting
% circles, given as a <matlab:doc('circleRegion') |circleRegion|> object.
% * A vector of real scalars specifying circulation strength around each
% circular boundary.
% * A vector of complex values specifing the locations of point vortices in
% the cirle domain.
% * A vector of real scalars giving the circulation strengths of the supplied
% point vortices.
% * Uniform background flow strength, a real scalar.
% * Uniform background flow angle, a real angle.
%
% Arguments to |flowRegion| can be given in a positional form. Any argument
% may be given as an empty array, with the exception of the circular
% islands, which must be given as an empty |circleRegion| object.

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

% Empty island argument example:
emptyIslands = circleRegion;


%%
% Arguments to |flowRegion| may also be given as name/value pairs. (See
% <matlab:doc('flowRegion') |flowRegion|> documentation for the list of
% property name arguments.)

pd2 = flowRegion(...
    'islands', circleRegion(circle(0, 1)), ...
    'islandCirculation', 0, ...
    'uniformStrength', uniformFlow, ...
    'uniformAngle', uniformAngle);


%%
% To facilitate evaluation of the potential function in the flow domain,
% one may use the class method |flowSamplePoints|. This computes a
% square grid of
% points in the domain with |npts| on a side which are fitted to a
% rectangular area |xylim|. Values for |xylim| are in the
% |[xmin, xmax, ymin, ymax]| format used with the <matlab:doc('axis')
% |axis|> command.
% The returned grid is a complex matrix of location values generated via
% the <matlab:doc('meshgrid') |meshgrid|> function
% with |NaN| values where there
% are circular obstacles or where points are within |vpad| of a point
% singularity. Default
% |vpad| value is |0.1|.
%
% In the following example we use the |plotbox| method built into
% |flowRegion| which returns suitible axis limits for plotting.

npts = 50;
xylim = plotbox(pd);
vpad = 0.3;
pts = flowSamplePoints(pd, npts, xylim, vpad);

clf
plot(pd)
hold on
plot(pts, 'b.')


%%
% Because properties of |flowRegion| may be updated individually, a check
% on the validity of the entire domain is only done on request. The way to
% do this is to use the |sanityCheck| method. This method returns an empty
% cell array if the domain is valid, otherwise the cell array contains
% strings describing what is wrong with the domain. This method is called
% automatically just before the |flowRegion| class constructor exits.

if ~isempty(sanityCheck(pd2))
    fprintf('Domain is not ok.\n')
end
pd2.islandCirculation = [1, 1];
chk = sanityCheck(pd2);
for k = 1:numel(chk)
    fprintf('%s\n', chk{k})
end
