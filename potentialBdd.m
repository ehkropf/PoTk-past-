classdef potentialBdd < complexPotential
%potentialBdd is potential in the bounded domain.

% E. Kropf, 2015

properties(SetAccess=protected)
    beta                    % Dipole for uniform flow.
    g0funs = {}             % Vortex functions.
    vjfuns = {}             % First kind integrals.
end

methods
    function W = potentialBdd(theDomain, varargin)
        W = W@complexPotential(varargin{:});
        if ~nargin
            return
        end
        
        if ~isa(theDomain, 'flowRegion')
            error('PoTk:InvalidArgument', 'Expected a flowRegion object.')
        end
        W.theDomain = theDomain;
        
        W = constructPotential(W);
    end
end

methods(Access=protected)
    function w = calcBdryCirc(W, zeta)
        % Circulation on the boundaries.
        
        w = complex(zeros(size(zeta)));
        
        gammav = W.theDomain.islandCirculation(2:end);
        for j = find(gammav')
            w = w + gammav(j)*W.vjfuns{j}(zeta);
        end
    end
    
    function w = calcPotential(W, zeta)
        % Combine potential flows.
        
        w = calcUniformFlow(W, zeta) ...
            + calcVortexCirc(W, zeta) ...
            + calcBdryCirc(W, zeta);
    end
    
    function w = calcUniformFlow(W, zeta)
        % Uniform flow potential.
        
        w = complex(zeros(size(zeta)));
        if isempty(W.beta) || W.theDomain.uniformStrength == 0
            return
        end
        
        % Do something here. Superclass method?
    end
    
    function w = calcVortexCirc(W, zeta)
        % Point vortex circulation.
        
        w = complex(zeros(size(zeta)));
        
        Gammav = W.theDomain.vortexCirculation;
        for k = find(Gammav')
            w = w + Gammav(k)*W.g0funs{k}(zeta);
        end
    end
    
    function W = constructPotential(W)
        % Call all the setup methods. In the proper order.
        
        W = setupBdryCirc(W);
        W = setupVortexCirc(W);
        W = setupUniformFlow(W);
    end

    function W = setupBdryCirc(W)
        % Construct the first kind integrals for boundary circulation
        % computation.
        % Construct them all regardless; then use to construct prime
        % functions.
        
        m = W.theDomain.m;
        D = W.theDomain.islands;
        
        W.vjfuns{1} = vjFirstKind(1, D);
        for j = 2:m
            W.vjfuns{j} = vjFirstKind(j, W.vjfuns{1});
        end
    end
    
    function W = setupUniformFlow(W)
        % Uniform flow depends on existence of point beta.
        
        if isempty(W.beta) || W.theDomain.uniformStrength == 0
            return
        end
        
        % Do something here? Superclass method?
    end
    
    function W = setupVortexCirc(W)
        % Construct Green's functions for vortices.
        % Extra boundary circulation goes to C0.
        
        alphav = W.theDomain.vortexLocation;
        n = numel(alphav);
        if n == 0
            return
        end
        
        for k = find(W.theDomain.vortexCirculation')
            W.g0funs{k} = greensC0(alphav(k), W.vjfuns{1});
        end
    end
end

end % classdef
