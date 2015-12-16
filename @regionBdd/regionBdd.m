classdef regionBdd < baseRegion
%regionBdd is the punctured unit disk region.

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
    
    function [zg, Lzg, axlim] = rectGrid(R, res, vpad)
        %[zg, Lzg, axlim] = rectGrid(R, res, vpad)
        %  Square point grid.
        
        if nargin < 3
            vpad = 0.1;
        end
        if nargin < 2
            res = 200;
        end
        
        [X, Y] = meshgrid(linspace(-1, 1, res));
        zg = complex(X, Y);
        
        cv = R.centers;
        rv = R.radii;
        Lzg = true(size(zg));
        Lzg(abs(zg) >= 1-eps(2)) = false;
        for j = 1:R.m
            Lzg(abs(zg - cv(j)) <= rv(j)+eps(2)) = false;
        end
        alphav = R.singularities;
        if vpad > 0
            for k = 1:numel(alphav)
                Lzg(abs(zg - alphav(k)) <= vpad) = false;
            end
        end
        
        if nargout > 2
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
