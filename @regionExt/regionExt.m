classdef regionExt < baseRegion
%regionExt is the region exterior to multiple circles.
%
% R = regionExt(dv, qv)
%   Constructs basic region defined by circle center vector dv and radius
%   vector qv.
% R = regionExt(dv, qv, circulation, singularities, singStrength, ...
%               uniformStrength, uniformAngle, dipole)
%   Constructs region by setting given parameters. The circulation vector
%   should have the same number of elements as circle centers given.
%   The uniform field strength and angle are real, finite values. The dipole
%   is the location of the uniform field dipole. Any parameter value may be
%   empty per MATLAB function convention.
% R = regionExt(dv, qv, 'parameterName', parameterValue, ...)
%   This is the name/value pair version of the constructor.

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

methods
    function R = regionExt(dv, qv, varargin)
        if ~nargin
            args = {};
        else
            args = [dv, qv, varargin];
        end
        R = R@baseRegion(args{:});
        if ~nargin
            return
        end
        
        sanityCheck(R)
    end
    
    function C = circleRegion(R)
        %C = circleRegion(R)
        %  Convert PoTk region to a CMT circle region.
        
        C = circleRegion(R.centers, R.radii);
    end
    
    function numberBoundaries(R)
        for j = 1:R.m+1
            c = R.centers(j);
            text(real(c), imag(c), num2str(j-1))
        end
    end
        
    function axlim = plotbox(R, scale)
        %axlim = plotbox(R, scale)
        
        if nargin < 2
            scale = [];
        end
        
        n = numel(R.singularities);
        if n && R.m
            axlim = cmt.plotbox([R.singularities; ...
                cmt.bb2z(boundbox(circleRegion(R)))], scale);
        elseif n
            axlim = cmt.plotbox(R.singularities, scale);
        elseif m
            axlim = plotbox(circleRegion(R), scale);
        else
            axlim = [-1, 1, -1, 1];
        end
    end
    
    function [zg, axlim] = rectGrid(R, res, axlim, vpad)
        %[zg, axlim] = rectGrid(R, res, axlim, vpad)
        %  Rectangular point grid.
        %  Argument axlim can be a valid axis array or a scalar where 1.0
        %  represents the smallest square to containt the circles and the
        %  singularity points.
        
        if nargin < 2
            res = 200;
        end
        if nargin < 3 || isempty(axlim)
            axlim = plotbox(R, 1.2);
        end
        if nargin < 4
            vpad = 0;
%             baselen = min(diff(axlim(1:2)), diff(axlim(3:4)));
%             vpad = 0.02*baselen;
        end
        
        if numel(axlim) == 1
            axlim = plotbox(R, axlim);
        end
        x = linspace(axlim(1), axlim(2), res);
        y = linspace(axlim(3), axlim(4), res);
        [X, Y] = meshgrid(x, y);
        zg = complex(X, Y);
        
        cv = R.centers;
        rv = R.radii;
        for j = 1:R.m+1
            zg(abs(zg - cv(j)) <= rv(j)+eps(2)) = nan;
        end
        alphav = R.singularities;
        if vpad > 0
            for k = 1:numel(alphav)
                zg(abs(zg - alphav(k)) <= vpad) = nan;
            end
        end
    end
    
    function [Rb, zeta] = regionBdd(R)
        %Construct bounded domain from exterior domain.
        %
        % [Rb, zeta] = regionBdd(R)
        % Constructs bounded region R from exterior region Re where zeta is
        % the Mobius transform by which this is done.
        
        if ~isa(R, 'regionExt')
            error(PoTk.ErrorTypeString.InvalidArgument, ...
                'Expected an exterior region object.')
        end
        if R.m < 1
            error(PoTk.ErrorTypeString.InvalidArgument, ...
                'Domain must have at least one boundary.')
        end
        
        % Mobius transformation.
        zeta = mobius(0, R.radii(1), 1, -R.centers(1));
        Rtmp = zeta(circleRegion(R));
        Rb = regionBdd(Rtmp.centers(2:end), Rtmp.radii(2:end));
        
        % Copy properties.
        Rb.circulation = R.circulation(2:end);
        Rb.singularities = zeta(R.singularities);
        Rb.dipole = zeta(R.dipole);        
        props = properties(R);
        props = props(5:end-1);
        for k = 1:numel(props)
            Rb.(props{k}) = R.(props{k});
        end
    end
end

methods(Access=protected)
    function m = mGetter(R)
        m = numel(R.centers) - 1;
    end
    
    function subSanityCheck(~)
        % Add checks to base.
        
        % Nothing yet.
    end
end

end
