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
    function W = setupCirculation(W)
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
    
    function W = setupVortices(W)
        % Construct Green's functions for vortices.
        % Extra boundary circulation goes to C0.
        
        alphav = W.theDomain.vortexLocation;
        Gammav = W.theDomain.vortexCirculation;
        n = numel(alphav);
        if n == 0
            return
        end
        
        for k = find(Gammav')
            W.g0funs{k} = greensC0(alphav(k), W.vjfuns{1});
        end
    end
end

end % classdef
