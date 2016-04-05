classdef potentialBdd < complexPotential
%potentialBdd is the bounded domain complex potential.
%
% W = potentialBdd(aDomain)
%   aDomain = regionBdd object.
%
% W = potentialBdd(aDomain, 'prop1', 'value1', ...)
%   Takes name/value pairs to set behaviour of the potential object. These
%   may be:
%   'useWaitBar' -- show a waitbar while doing calculations (default = false).
%   'streamWithField' -- plot stream lines with velocity field
%     (default = true).
%
% See also regionBdd.

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

properties(SetAccess=protected)
    g0funs = {}             % Vortex functions.
    vjfuns = {}             % First kind integrals.
    dG0db = {}              % Green's derivatives wrt parameter.
end

properties(Access=protected)
    defaultPlotScale = 1
end

methods
    function W = potentialBdd(theDomain, varargin)
        subbar = [];
        if numel(varargin) > 0 ...
                && isa(varargin{1}, 'PoG.subbar')
            subbar = varargin{1};
            varargin = varargin(2:end);
        end
        
        W = W@complexPotential(varargin{:});
        if ~nargin
            return
        end
        
        if ~isa(theDomain, 'regionBdd')
            error(PoTk.ErrorTypeString.InvalidArgument, ...
                'Expected a "regionBdd" object.')
        end
        W.theDomain = theDomain;
        
        if ~isempty(subbar)
            W.useWaitBar = true;
            W.waitBar = subbar;
        end
        
        W = constructPotential(W);
    end
end

methods(Access=protected)
    %%%%% Calculation.
    function w = calcPotential(W, zeta, subbar)
        % Combine potential flows.
        
        if nargin > 2 && isa(subbar, 'PoG.subbar')
            W.waitBar = subbar;
        end
        
        W = waitbarInitialize(W, 'Computing potential values');
        w = calcBdryCirc(W, zeta);
        w = w + calcVortexCirc(W, zeta);
        w = w + calcUniformFlow(W, zeta);
        waitbarRelease(W);
    end
    
    function w = calcBdryCirc(W, zeta)
        % Circulation on the boundaries.
        
        w = complex(zeros(size(zeta)));
        
        gammav = W.theDomain.circulation;
        m = numel(gammav);
        for j = find(gammav')
            waitbarUpdate(W, (j-1)/m/3, ...
                sprintf('Boundary contribution (%d/%d) at %d points', ...
                j, m, numel(zeta)))
            w = w + gammav(j)*W.vjfuns{j}(zeta);
        end
    end
    
    function w = calcVortexCirc(W, zeta)
        % Point vortex circulation.
        
        w = complex(zeros(size(zeta)));
        
        Gammav = W.theDomain.singStrength;
        n = numel(Gammav);
        for k = find(Gammav')
            waitbarUpdate(W, (1 + (k-1)/n)/3, ...
                sprintf('Singularity contribution (%d/%d) at %d points', ...
                k, n, numel(zeta)))
            w = w + Gammav(k)*W.g0funs{k}(zeta);
        end
    end
    
    function w = calcUniformFlow(W, zeta)
        % Uniform flow potential.
        
        waitbarUpdate(W, 2/3, 'Uniform field contribution')
        
        beta = W.theDomain.dipole;
        U = W.theDomain.uniformStrength;
        Chi = -W.theDomain.uniformAngle;
        if isempty(beta) || U == 0
            w = complex(zeros(size(zeta)));
            return
        end
        
        switch W.theDomain.m
            case 0
                w = U*zeta*exp(-1i*Chi);
                return
                
            case 1
                zMinusBeta = zeta - beta;
                w = U*(exp(-1i*Chi)./zMinusBeta ...
                    + exp(1i*Chi)*zMinusBeta);
                
            otherwise
                w = complex(zeros(size(zeta)));
                h = W.hCenterDiff;
                if mod(Chi, pi) > eps(pi)
                    % There is a horizontal component.
                    w = w + ...
                        (W.dG0db{1}(zeta) - W.dG0db{2}(zeta))/h*sin(Chi);
                end
                if mod(Chi + pi/2, pi) > eps(pi)
                    % There is a vertical component.
                    w = w + ...
                        (W.dG0db{3}(zeta) - W.dG0db{4}(zeta))/h*cos(Chi);
                end
                w = -4*pi*U*w;
        end
    end
    
    %%%%% Construction.
    function W = constructPotential(W)
        % Call all the setup methods. In the proper order.
        
        W = waitbarInitialize(W, 'Configuring potential functions');
        W = setupBdryCirc(W);
        W = setupVortexCirc(W);
        W = setupUniformFlow(W);
        W = waitbarRelease(W);
    end

    function W = setupBdryCirc(W)
        % Construct the first kind integrals for boundary circulation
        % computation.
        % Construct them all regardless; then use to construct prime
        % functions.
        
        m = W.theDomain.m;
        dv = W.theDomain.centers;
        qv = W.theDomain.radii;
        
        waitbarUpdate(W, 0, sprintf('Boundary part (1/%d)', m))
        W.vjfuns{1} = vjFirstKind(1, skpDomain(dv, qv));
        for j = 2:m
            waitbarUpdate(W, (j-1)/m/3, ...
                sprintf('Boundary part (%d/%d)', j, m))
            W.vjfuns{j} = vjFirstKind(j, W.vjfuns{1});
        end
    end
    
    function W = setupVortexCirc(W)
        % Construct Green's functions for vortices.
        % Extra boundary circulation goes to C0.
        
        alphav = W.theDomain.singularities;
        n = numel(alphav);
        if n == 0
            return
        end
        
        for k = find(W.theDomain.singStrength')
            waitbarUpdate(W, (1 + (k-1)/n)/3, ...
                sprintf('Singularity part (%d/%d)', k, n))
            W.g0funs{k} = greensC0(alphav(k), W.vjfuns{1});
        end
    end
    
    function W = setupUniformFlow(W)
        % Uniform flow depends on existence of point beta.
        
        waitbarUpdate(W, 2/3, 'Unform field part')
        
        beta = W.theDomain.dipole;
        if isempty(beta) || W.theDomain.m < 2 ...
                || W.theDomain.uniformStrength == 0
            return
        end
        
        Chi = W.theDomain.uniformAngle;
        h = W.hCenterDiff;
        db = beta + 0.5*h*[1, -1, 1i, -1i];
        W.dG0db = cell(4, 1);
        
        if mod(Chi, pi) > eps(pi)
            % There is a horizontal component.
            W.dG0db{1} = greensC0(db(1), W.vjfuns{1});
            W.dG0db{2} = greensC0(db(2), W.vjfuns{1});
        end
        if mod(Chi + pi/2, pi) > eps(pi)
            % There is a vertical component.
            W.dG0db{3} = greensC0(db(3), W.vjfuns{1});
            W.dG0db{4} = greensC0(db(4), W.vjfuns{1});
        end
    end
end

end % classdef
