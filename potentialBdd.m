classdef potentialBdd < complexPotential
%potentialBdd is potential in the bounded domain.

% E. Kropf, 2015

properties(SetAccess=protected)
    beta                    % Dipole for uniform flow.
end

methods
    function W = potentialBdd(theDomain, varargin)
        if ~nargin
            return
        end
        W = W@complexPotential(varargin{:});
        
        if ~isa(theDomain, 'flowRegion')
            error('PoTk:InvalidArgument', 'Expected a flowRegion object.')
        end
        W.primeDomain = theDomain;
        
        W = constructPotential(W);
    end
end

methods(Access=protected)
    % Computation needs only G_0 at the vortices and v_j for boundary
    % circulation. No need for intermediate "beta" point to "transfer"
    % circulation.
    
    % Then setup only needs to construct G_0 and get the v_j.
end

end % classdef
