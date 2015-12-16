classdef regionExt < baseRegion
%regionExt is the region exterior to multiple circles.

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
    
    function [zg, Lzg, axlim] = rectGrid(R, res, axlim, vpad)
        %[zg, Lzg, axlim] = rectGrid(R, res, axlim, vpad)
        %  Rectangular point grid.
        
        if nargin < 4
            vpad = 0.1;
        end
        if nargin < 3 || isempty(axlim)
            axlim = plotbox(R, 1.2);
        end
        if nargin < 2
            res = 200;
        end
        
        x = linspace(axlim(1), axlim(2), res);
        y = linspace(axlim(3), axlim(4), res);
        [X, Y] = meshgrid(x, y);
        zg = complex(X, Y);
        
        cv = R.centers;
        rv = R.radii;
        Lzg = true(size(zg));
        for j = 1:R.m
            Lzg(abs(zg - cv(j)) <= rv(j)+eps(2)) = false;
        end
        alphav = R.singularities;
        if vpad > 0
            for k = 1:numel(alphav)
                Lzg(abs(zg - alphav(k)) <= vpad) = false;
            end
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
