classdef(Abstract) complexPotential
%complexPotential is the abstract base class for complex potentials.
%
%  This abstract class may not be instantiated directly.
%
%  Constructor:
%  P = complexPotential('name1', value1, ...)
%  Provide name/value pair processing for properties
%    useWaitBar
%    streamWithField


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
    useWaitBar = false      % Show wait bar during calculations.
    streamWithField = true  % Plot stream lines with field lines.
end

properties(SetAccess=protected)
    theDomain               % (r/o) basic domain.
end

properties(Abstract,Access=protected)
    defaultPlotScale        % Default scale for plotboxes.
end

properties(Access=protected)
    zetaFun                 % (p) Map from input to prime domain.

    hCenterDiff = 1e-3      % (p) Centred difference h value.
    numSamplePts = 75       % (p) sqrt(#) of sample points for stream lines in
                            % one dimension.
    numPlotPts = 200        % (p) sqrt(#) of points to use for plotting 
                            % stream lines.
    numVectorPts = 20       % (p) sqrt(#) of points for vortex refinement.
    waitBar                 % Waitbar object.    
end

properties(SetAccess=private,GetAccess=protected)
    varPropList = {'useWaitBar', 'streamWithField'}
end

methods
    % Constructor is protected method; see below.
    % FIXME: If this class is abstract, a protected constructor doesn't
    % make too much sense.
    
    function disp(W)
        cname = class(W);
        fprintf(['  <a href="matlab:helpPopup %s">%s</a> ' ...
            'object with properties:\n'], cname, cname)
        fprintf('\n   input domain has %d boundary components\n', W.theDomain.m)
%         fprintf('   number of stored prime functions: %d\n', ...
%             numel(W.SKav) + numel(W.SKbv) + numel(W.SKbxy))
        fprintf('   show wait bar while computing: %s\n', ...
            PoTk.logical2str(W.useWaitBar))
        fprintf('   plot stream lines with velocity field: %s\n', ...
            PoTk.logical2str(W.streamWithField))
        fprintf('\n')
    end
    
    function w = feval(W, z)
        % Provide feval capability.
        
        w = calcPotential(W, z);
    end
    
    function W = set.useWaitBar(W, value)
        W.useWaitBar = validateLogical(W, value);
    end
    
    function W = set.streamWithField(W, value)
        W.streamWithField = validateLogical(W, value);
    end
    
    function plot(W, varargin)
        % Auto plot potential information; default is to plot stream lines.
        %
        % plot(W)
        % plot(W, 'streamLines', ...)
        %   Plots the streamlines for the potential along with the domain.
        %   This calls plotStreamLines(W, ...); extra arguments are passed to
        %   this function.
        %   W = potential object.
        %
        % plot(W, 'velocityField', ...)
        %   Plots the velocity vector field for the potential along with
        %   the domain. This calls plotVelocityField(W, ...); extra
        %   arguments are passed to this function.
        
        validStrings = {'streamLines', 'velocityField'};
        defaultPlot = validStrings{1};
        if nargin > 1
            try
                vstr = validatestring(varargin{1}, validStrings);
                varargin = varargin(2:end);
            catch
                vstr = defaultPlot;
            end
        else
            vstr = defaultPlot;
        end
        
        switch vstr
            case validStrings{1}
                plotStreamLines(W, varargin{:})
            case validStrings{2}
                plotVelocityField(W, varargin{:})
        end
    end
     
    function plotStreamLines(W, varargin)
        % Plot stream lines along with domain.
        % Plots to current axes, or chooses axes from physical domain.
        % Extra arguments are passed to the MATLAB built in contour.
        %
        % plotStreamLines(W, ...)
        %   W = potential object.
        
        washold = ishold;
        if ~washold
            cah = newplot;
            axis off
            cmtplot.whitefigure(cah)
            axlim = plotbox(W.theDomain, W.defaultPlotScale);
            hold on
        else
            axlim = axis;
        end
        axscale = max(diff(axlim(1:2)), diff(axlim(3:4)));
        
        % Calculate.        
        zs = rectGrid(W.theDomain, ...
            W.numPlotPts, axlim, 0.005*axscale);
        Wz = feval(W, zs);
        
        % Ignore warning if no flow (F is constant).
        warning('off', 'MATLAB:contour:ConstantData')
        contour(real(zs), imag(zs), imag(Wz), ...
            20, 'lineColor', [137, 196, 232]/255, varargin{:})
        warning('on', 'MATLAB:contour:ConstantData')

        % Draw.
        hold on        
        plot(W.theDomain)
        
        if ~washold
            numberBoundaries(W.theDomain)
            axis(axlim)
            aspectequal
            hold off
        end
    end
    
    function plotVelocityField(W, varargin)
        % Plot a vector field from the potential along with the domain.
        % Plots to current axes, or chooses axes from physical domain.
        % Extra arguments are passed to the MATLAB built in quiver.
        %
        % plotVelocityField(W, ...)
        %   W = potential object.
        
        cah = newplot;
        washold = ishold;
        if ~washold
            cmtplot.whitefigure(cah)
            axlim = plotbox(W.theDomain, 1.5);
            axis(axlim)
        else
            axlim = axis;
        end
        
        xpad = 0.01*diff(axlim([1,2]));
        ypad = 0.01*diff(axlim([3,4]));
        vlim = [axlim(1)+xpad, axlim(2)-xpad, axlim(3)+ypad, axlim(4)-ypad];
        vscale = max(diff(vlim(1:2)), diff(vlim(3:4)));
        zs = flowSamplePoints(W.theDomain, W.numVectorPts, ...
            vlim, 0.03*vscale);
        zs = zs(~isnan(zs));
        V = velocityField(W, zs);
        
        hold on
        doStreamLines = true;
        if nargin > 1 && strcmp(varargin{1}, 'noStreamLines')
            doStreamLines = false;
            varargin = varargin(2:end);
        end
        if W.streamWithField && doStreamLines
            plot(W)
        else
            plot(W.theDomain)
        end
        quiver(real(zs), imag(zs), real(V), imag(V), ...
            'color', [0.078, 0.169, 0.549], varargin{:})
        
        if ~washold
            numberIslandsAndVortices(W.theDomain)
            aspectequal
            axis off
            hold off
        end
    end
   
    function out = subsref(W, S)
        % Provide w = W(zeta) functionality.
        % Calls w = feval(W, zeta).
        
        if numel(S) == 1 && strcmp(S.type, '()')
            out = calcPotential(W, S.subs{:});
        else
            out = builtin('subsref', W, S);
        end
    end
    
    function V = velocityField(W, Z)
        % Calculate the velocity vector field values at points Z.
        % Uses the central difference method for taking derivatives of the
        % potential. It's assumed the sample
        % points are more that 0.5*abs(h) away from islands or vortices.
        %
        % V = velocityField(W, Z)
        %   W = potential object.
        %   Z = Points at which to find the velocity.
        %
        %   V = Calculated velocity at Z.
        
        h = W.hCenterDiff;
        
        wasWaitBar = W.useWaitBar;        
        if wasWaitBar
            W.useWaitBar = false;
            msg = sprintf(...
                ['Calculating derivatives of the potential at ' ...
                '%d points. (%%d/4)'], numel(Z));
            wbh = waitbar(0, sprintf(msg, 1), ...
                'name', 'Calculate flow field.');
        end
        u = imag(feval(W, Z + 0.5i*h))/h;
        
        if wasWaitBar
            waitbar(1/4, wbh, sprintf(msg, 2))
        end
        u = u - imag(feval(W, Z - 0.5i*h))/h;
        
        if wasWaitBar
            waitbar(1/2, wbh, sprintf(msg, 3))
        end
        v = -imag(feval(W, Z + 0.5*h))/h;
        
        if wasWaitBar
            waitbar(3/4, wbh, sprintf(msg, 4))
        end
        v = v + imag(feval(W, Z - 0.5*h))/h;
        
        if wasWaitBar
            delete(wbh)
            drawnow
        end
        
        V = u + 1i*v;
    end    
end

methods(Hidden)
    function [w, z] = computePotentialField(W, numPts, xyLim)
        % Compute potential field for stream line plot.

        warning('PoTk:deprecatedFunction', ...
            ['The use of this function is deprecated. ' ...
            'It may be removed without warning in future versions.'])
        
        ns = W.numSamplePts;
        if numPts < ns
            ns = numPts;
        end

        if ns ~= numPts
            % Basic sample grid.
            zs = flowSamplePoints(W.theDomain, ns, xyLim);
        
            % Finer grid around vortices.
            vortices = W.theDomain.vortexLocation;
            n = numel(vortices);
            if n > 1
                np = 16;
                [X, Y] = meshgrid(linspace(-1, 1, np));
                zp(1:n*size(X, 1), 1:size(X, 2)) = 1i;
                for k = 1:n
                    zp((k-1)*np+(1:np),:) = vortices(k) + complex(X, Y)*0.3;
                end
            else
                zp = [];
            end
        
            % Evaluate.
            zz = [zs(:); zp(:)];
            w = complex(nan(size(zz)));
            L = ~isnan(zz);
            w(L) = feval(W, zz(L));
        
            S = scatteredInterpolant(...
                real(zz(L)), imag(zz(L)), w(L), 'natural');
            z = flowSamplePoints(W.theDomain, numPts, xyLim, 0.05);
            L = ~isnan(z);
            w = complex(nan(size(z)));
            w(L) = S(real(z(L)), imag(z(L)));
        else
            z = flowSamplePoints(W.theDomain, numPts, xyLim, 0);
            w = complex(nan(size(z)));
            L = ~isnan(z);
            w(L) = feval(W, z(L));
        end
        w(~L) = max(w(L));
    end
    
    function [V, z] = computeVelocityField(W, numPts, xyLim)
        % Compute velocity field for velocity plot.
        
        warning('PoTk:deprecatedFunction', ...
            ['The use of this function is deprecated. ' ...
            'It may be removed without warning in future versions.'])
        
        ns = W.numVectorPts;
        if numPts < ns
            ns = numPts;
        end
        
        % Basic grid plane.
        xyscale = max(diff(xyLim(1:2)), diff(xyLim(3:4)));
        z = flowSamplePoints(W.theDomain, ns, xyLim, 0.03*xyscale);
                
        % Compute.
        z = z(~isnan(z));
        V = velocityField(W, z);
    end
end

methods(Access=protected)
    function cpot = complexPotential(varargin)
        % Provide name/value pair processing for properties
        %   useWaitBar
        %   streamWithField
        
        if ~nargin
            return
        end
                
        n = numel(varargin);
        for j = 1:2:n-1
            if ~ischar(varargin{j})
                error('PoTk:InvalidArgument', ...
                    'Expected name/value pairs following physical domain.')
            end
            try
                vstr = validatestring(varargin{j}, cpot.varPropList);
            catch me
                error(struct(...
                    'identifier', 'PoTk:invalidArgument', ...
                    'message', ...
                    sprintf('Property "%s" not recognized.', varargin{j}), ...
                    'stack', me.stack))
            end
            cpot.(vstr) = varargin{j+1};
        end
    end
        
    function tf = validateLogical(~, value)
        % Logical value validation for set functions.
        
        tf = logical(value);
        if ~islogical(tf)
            error('PoTk:InvalidArgument', 'Excpected logical value.')
        end
    end
    
    function W = waitbarInitialize(W, name, msg)
        if ~W.useWaitBar
            return
        end
        if nargin < 3
            msg = '';
        end
        
        if ~isempty(W.waitBar) && isa(W.waitBar, 'PoG.subbar')
            update(W.waitBar, 0, msg);
        else
            W.waitBar = PoG.waitbar(name, msg);
        end
    end
    
    function waitbarUpdate(W, x, msg)
        if ~W.useWaitBar || isempty(W.waitBar)
            return
        end
        update(W.waitBar, x, msg)
    end
    
    function W = waitbarRelease(W)
        if ~W.useWaitBar || isempty(W.waitBar)
            return
        end
        release(W.waitBar)
        W.waitBar = [];
        drawnow
    end
    
    %%%%% Abstract protected methods.
    W = constructPotential(W)
    w = calcPotential(W, zeta)
end

methods(Access=protected,Static)
    function a = annulus(s, r)
        % Calculate the value a so that the circle with center s and radius
        % r (with |s| < 1 and 0 < |s| + r < 1) is mapped to the interior
        % annulus boundary, given the outer unit disk boundary applies to
        % both circles. The value a is used in the map
        %
        %           z - a
        % A(z) = -------------
        %        1 - z*conj(a)
        
        % Put circle on positive real axis
        k = exp(1i*angle(s));
        s2 = abs(s)^2;
        r2 = r^2;        
        
        a = (1 + s2 - r2 - ...
            sqrt((s2 - r2)^2 - 2*(s2 + r2) + 1))/(2*sqrt(s2));
        assert(a >= 0, 'Unexpected non-positive value detected.')
        a = k*a;
    end
end

end
