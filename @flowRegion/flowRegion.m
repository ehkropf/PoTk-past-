classdef flowRegion
%flowRegion represents a flow domain for the PoTk.
%
% fd = flowRegion(islands, islandCirculation, ...
%                 vortexLocation, vortexCirculation, ...
%                 uniformStrength, unformAngle)
% Constructs a flow domain object. All values are optional. See below for
% the input descriptions.
%
% fd = flowRegion('property1', value1, 'property2', value2, ...)
% Same as above using the name/value pair argument scheme.
%
% Inputs:
% (The following can be used in the name/value pair scheme.)
%   islands - a circleRegion object.
%   islandCirculation - a vector of real scalars, must be same length as
%       islands.
%   vortexLocation - a vector of complex values.
%   vortexCirculation - a vector of real scalars the same length as
%       vortexCirculation.
%   uniformStrength - a real scalar value.
%   uniformAngle - a real scalar value which will be nomalized to [0,2*pi).
%
% See also circleRegion.

% Everett Kropf, 2015
% 
% This file is part of the Potential Toolkit (PoTk).
% 
% PoTk is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% PoTk is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with PoTk.  If not, see <http://www.gnu.org/licenses/>.

properties
    islands = circleRegion      % Circular obstructions.
    islandCirculation           % Flow circulation on island boundary.
    vortexLocation              % Point vortex location vector.
    vortexCirculation           % Vortex circulation strength vector.
    uniformStrength = 0         % Unform background flow strength.
    uniformAngle = 0            % Unform flow angle.
    dipole = complex(inf)       % Uniform flow dipole location.
end

properties(Dependent)
    m                           % (Read-only) Domain connectivity.
end

methods
    function fd = flowRegion(varargin)
        if ~nargin
            return
        end
        
        props = properties(flowRegion);
        if isa(varargin{1}, 'circleRegion')
            for j = 1:nargin
                fd.(props{j}) = varargin{j};
            end
        elseif ischar(varargin{1})
            for j = 1:2:nargin-1
                try
                    vstr = validatestring(varargin{j}, props);
                    fd.(vstr) = varargin{j+1};
                catch err
                    switch err.identifier
                        case 'MATLAB:unrecognizedStringChoice'
                            errmsg = sprintf(...
                                'Parameter ''%s'' unrecognized.', varargin{j});
                        case 'MATLAB:badsubsript'
                            errmsg = sprintf(...
                                'Value expected for parameter ''%s''.', ...
                                varargin{j});
                        otherwise
                            rethrow(err)
                    end
                    error('PoTk:InvalidArgument', '%s', errmsg)
                end
            end
        else
            %error
        end
        
        msg = sanityCheck(fd);
        if ~isempty(msg)
            msg = sprintf('%s', msg{:});
            error('PoTk:UndefinedState', ...
                'flowRegion is in an undefined state:\n%s', msg)
        end
    end

    function disp(fd)
        fprintf(['  <a href="matlab:helpPopup flowRegion">' ...
            'flowRegion</a> object with properties:\n'])
        fprintf('\n  connectivity = %d\n', fd.m)
        str = strtrim(evalc('disp(fd.islandCirculation'')'));
        fprintf('  island circulation: %s\n', str)
        fprintf('  vortex count = %d\n', numel(fd.vortexLocation))
        str = strtrim(evalc('disp(fd.vortexCirculation'')'));
        fprintf('  vortex circulation: %s\n', str)
        if fd.uniformStrength
            str = strtrim(evalc('disp(fd.uniformStrength)'));
            fprintf('  background flow strength = %s\n', str)
            str = strtrim(evalc('disp(fd.uniformAngle/pi)'));
            fprintf('  background flow angle = %s*pi\n', str)
        else
            fprintf('  no background flow\n')
        end
        fprintf('\n')
    end
    
    function Z = flowSamplePoints(fd, npts, xyLim, vpad)
        % Return a rectangular grid of points in the domain.
        % Grid points inside circle objects are replaced with NaN values.
        % Grid points within a radius of vpad of a vortex are also set to
        % NaN.
        %
        % Z = flowSamplePoints(fd, npts, xyLim, vpad)
        %   fd = flowRegion object.
        %   npts = number of grid points on a side to use.
        %   xyLim = rectangular area in which to fit grid, format is
        %     [xmin, xmax, ymin, ymax].
        %   vpad = (optional) size of radius around vortices
        %     (default = 0.1).
        %
        % Grid points are simply created with meshgrid and then replaced
        % with NaN where approprate.
        
        nZ = npts;
        if nargin < 4
            vpad = 0.1;
        end
        if nargin < 3 || isempty(xyLim)
            xyLim = plotbox(fd, 1.2);
        end
        
        x = linspace(xyLim(1), xyLim(2), nZ);
        y = linspace(xyLim(3), xyLim(4), nZ);
        [X, Y] = meshgrid(x, y);
        
        Z = complex(X, Y);
        cv = fd.islands.centers;
        rv = fd.islands.radii;
        for j = 1:numel(cv)
            Z(abs(Z - cv(j)) <= rv(j)) = nan;
        end
        alphas = fd.vortexLocation;
        if vpad > 0
            for k = 1:numel(alphas)
                Z(abs(Z - alphas(k)) <= vpad) = nan;
            end
        end
    end
    
    function numberIslandsAndVortices(fd)
        for j = 1:numel(fd.islands.circles)
            c = fd.islands.centers(j);
            text(real(c), imag(c), num2str(j))
        end
        d = axis;
        d = 0.01*max(diff(d(1:2)), diff(d(3:4)));
        for k = 1:numel(fd.vortexLocation)
            v = fd.vortexLocation(k);
            text(real(v) + d, imag(v) + d, num2str(k), 'color', 'm')
        end
    end
    
    function plot(fd)
        % Draw the flowRegion.
        %
        % plot(fd) where fd is a flowRegion object.
        
        cah = newplot;
        washold = ishold;
        
        hold on
        fill(inv(fd.islands))
        plot(fd.islands)
        z = fd.vortexLocation;
        plot(real(z), imag(z), 'k.', 'markersize', 18)
        
        if ~washold
            cmtplot.whitefigure(cah)
            axis(plotbox(fd))
            aspectequal
            axis off
            hold off
        end
    end

    function xylim = plotbox(fd, scale)
        % Return rectangular box around domain for plot axis.
        %
        % xylim = plotbox(fd, scale)
        %   fd = flowRegion object.
        %   scale = (optional) scalar specifying extra space around
        %     bounding box for the plotbox.
        %
        %   xylim = array suitible to use in the axis() command. Format
        %     is [xmin, xmax, ymin, ymax].
        
        if nargin < 2
            scale = [];
        end
        
        n = numel(fd.vortexLocation);
        m = fd.m;
        if n && m
            xylim = cmt.plotbox([fd.vortexLocation; ...
                cmt.bb2z(boundbox(fd.islands))], scale);
        elseif n
            xylim = cmt.plotbox(fd.vortexLocation, scale);
            if all(xylim == 0)
                xylim = [-1, 1, -1, 1];
            end
        elseif m
            xylim = plotbox(fd.islands, scale);
        else
            xylim = [-1, 1, -1, 1];
        end
    end
    
    function msg = sanityCheck(fd)
        % Check that flow domain makes sense as configured.
        % Domain must not have intersecting circles, vortices inside
        % circles, and must have circulation specified for every circle and
        % vortex.
        %
        % msg = sanityCheck(fd)
        %   fd = flowRegion object.
        %
        %   msg = cell array of strings describing problmes with domain.
        %     The array is empty if the domain passes the checks.
        
        msg = {};
        
        % # of islands and circluations must be the same.
        if ~isempty(fd.islands)
            n = numel(boundary(fd.islands));
        else
            n = 0;
        end
        if n ~= numel(fd.islandCirculation)
            msg{end+1} = sprintf(...
                'Number of islands and circulation values do not match.\n');
        end
        
        % # of vortices and circulations must be the same.
        if numel(fd.vortexLocation) ~= numel(fd.vortexCirculation)
            msg{end+1} = sprintf(...
                'Number of vortices and vortex strength values are not equal.\n');
        end
        
        % No vortices inside islands.
        if ~isempty(fd.islands) && ~isempty(fd.vortexLocation)
            alpha = fd.vortexLocation;
            for k = 1:numel(alpha)
                if any(abs(fd.islands.centers - alpha(k)) < fd.islands.radii)
                    msg{end+1} = sprintf(...
                        'Vortex detected inside island boundary.\n'); %#ok<AGROW>
                    break
                end
            end
        end
    end
    
    %----------------------------------------------------------------
    function m = get.m(fd)
        m = fd.islands.m;
    end
    
    function fd = set.islands(fd, cregion)
        if ~isa(cregion, 'circleRegion')
            errorInvalidArg(fd, ...
                'The ''islands'' property expects a circleRegion object.')
        end
        fd.islands = cregion;
    end
    
    function fd = set.islandCirculation(fd, circv)
        circv = circv(:);
        if any(~isRealAndFinite(fd, circv))
            errorInvalidArg(fd, ...
                'Island circulation values must be real and finite.')
        end
        fd.islandCirculation = circv;
    end
    
    function fd = set.vortexLocation(fd, locv)
        locv = locv(:);
        if any(~isFiniteNumber(fd, locv))
            errorInvalidArg(fd, ...
                'Vortex locations must be finite numbers')
        end
        fd.vortexLocation = locv;
    end
    
    function fd = set.vortexCirculation(fd, circv)
        circv = circv(:);
        if any(~isRealAndFinite(fd, circv))
            errorInvalidArg(fd, ...
                'Vortex circulation values must be real and finite.')
        end
        fd.vortexCirculation = circv;
    end
    
    function fd = set.uniformStrength(fd, stval)
        if numel(stval) > 1 || ~isRealAndFinite(fd, stval)
            errorInvalidArg(fd, ...
                'Uniform flow strength must be a real, finite, scalar number.')
        end
        fd.uniformStrength = stval;
    end
    
    function fd = set.uniformAngle(fd, arg)
        if numel(arg) > 1 || ~isRealAndFinite(fd, arg)
            errorInvalidArg(fd, ...
                'Uniform flow angle must be a real, finite, scalar number.')
        end
        fd.uniformAngle = arg;
    end
end

methods(Access=protected)
    function errorInvalidArg(~, errmsg, varargin)
        % Common invalid argument error generator.
        
        error('PoTk:InvalidArgument', errmsg, varargin{:})
    end
    
    function tf = isFiniteNumber(~, vec)
        % Numeric & finite check.
        
        tf = isnumeric(vec) & isfinite(vec);
    end
    
    function tf = isRealAndFinite(fd, vec)
        % Numeric & finite & real check.
        
        tf = isFiniteNumber(fd, vec) & imag(vec) == 0;
    end
end

end
