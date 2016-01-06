classdef regionBdd < baseRegion
%regionBdd is the punctured unit disk region.
%
% R = regionBdd(dv, qv)
%   Constructs basic region defined by circle center vector dv and radius
%   vector qv.
% R = regionBdd(dv, qv, circulation, singularities, singStrength, ...
%               uniformStrength, uniformAngle, dipole)
%   Constructs region by setting given parameters. The circulation vector
%   should have one more element at the start of the vector than the number
%   of circle centers (circulation on the unit circle). The singularities
%   and singStrength vectors must be the same size. The uniform field 
%   strength and angle are real, finite values. The dipole is the location 
%   of the uniform field dipole. Any parameter value may be empty per MATLAB
%   function convention.
% R = regionBdd(dv, qv, 'parameterName', parameterValue, ...)
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
    function R = regionBdd(dv, qv, varargin)
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
        
        C = circleRegion([0; R.centers], [1; R.radii]);
    end
    
    function numberBoundaries(R)
        for j = 1:R.m
            c = R.centers(j);
            text(real(c), imag(c), num2str(j))
        end
    end
    
    function D = skpDomain(R)
        % Convert to SKPrime domain object.
        
        D = skpDomain(R.centers, R.radii);
    end
    
    function xylim = plotbox(~, scale)
        if nargin < 2
            scale = 1;
        end
        xylim = scale*[-1, 1, -1, 1];
    end
    
    function [zg, axlim] = rectGrid(R, res, vpad)
        %[zg, axlim] = rectGrid(R, res, vpad)
        %  Square point grid.
        
        if nargin < 2
            res = 200;
        end
        if nargin < 3
            vpad = 0; %0.005;
        end
        
        [X, Y] = meshgrid(linspace(-1, 1, res));
        zg = complex(X, Y);
        
        dv = R.centers;
        qv = R.radii;
        zg(abs(zg) >= 1-eps(2)) = nan;
        for j = 1:R.m
            zg(abs(zg - dv(j)) <= qv(j)+eps(2)) = nan;
        end
        alphav = R.singularities;
        if vpad > 0
            for k = 1:numel(alphav)
                zg(abs(zg - alphav(k)) <= vpad) = nan;
            end
        end
        
        if nargout > 1
            axlim = [-1, 1, -1, 1];
        end
    end
end

methods(Access=protected)
    function m = mGetter(R)
        m = numel(R.centers);
    end
    
    function subSanityCheck(~)
        % Add checks to base set.
        
        % Nothing yet.
    end
end

        end
    end
end

end
