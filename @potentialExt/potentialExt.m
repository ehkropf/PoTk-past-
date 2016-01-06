classdef potentialExt < complexPotential
%potentialExt is the unbounded domain complex potential.
%
% W = potentialExt(aDomain)
%   aDomain = regionExt object.
%
% W = potentialExt(aDomain, 'prop1', 'value1', ...)
%   Takes name/value pairs to set behaviour of the potential object. These
%   may be:
%   'useWaitBar' -- show a waitbar while doing calculations (default = false).
%   'streamWithField' -- plot stream lines with velocity field
%     (default = true).
%
% See also regionExt.

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
    bddDomain               % (r/o) Bounded domain.
    bddPotential            % (r/o) Bounded potential.
    zetaf                   % (r/o) Mobius transform unbounded -> bounded.
    beta = inf              % (r/o) Image of infinity in bounded domain.
    greensBeta              % (r/o) Green's function wrt C0 at beta.
end

properties(Access=protected)
end

methods
    function W = potentialExt(theDomain, varargin)
        W = W@complexPotential(varargin{:});
        if ~nargin
            return
        end
        
        if ~isa(theDomain, 'regionExt')
            error(PoTk.ErrorTypeString.InvalidArgument, ...
                'Expected a "regionExt" object.')
        end        
        W.theDomain = theDomain;
        
        W = constructPotential(W);
    end
end

methods(Access=protected)
    %%%%% Calculation.
    function w = calcPotential(W, z)
        % Combine with internal calculation.
        
        zeta = W.zetaf(z);
        w = calcBounded(W, zeta);
        w = w - calcExternalPart(W, zeta);
    end
    
    function w = calcBounded(W, zeta)
        % Bounded potential part.
        
        w = W.bddPotential(zeta);
    end
    
    function w = calcExternalPart(W, zeta)
        % Calculate extra part for external flow.
        
        extSums = sum(W.theDomain.singStrength) ...
            + sum(W.theDomain.circulation);
        if extSums == 0
            w = complex(zeros(size(zeta)));
            return
        end
        w = extSums*W.greensBeta(zeta);
    end
    
    %%%%% Construction.
    function W = constructPotential(W)
        
        W = setupBoundedRegion(W);
        W = setupBoundedPotential(W);
        W = setupGreensFunction(W);
    end
    
    function W = setupBoundedRegion(W)
        % Mobius transform of exterior domain.
        
        if W.theDomain.m > 0
            [Db, zeta] = regionBdd(W.theDomain);
            
            W.beta = pole(inv(zeta));
            W.zetaf = zeta;
            W.bddDomain = Db;
        else
            % No equivalent bounded domain.
            W.zetaf = mobius(1, 0, 0, 1);
        end
    end
    
    function W = setupBoundedPotential(W)
        % Get the bounded potential.
        
        W.bddPotential = potentialBdd(W.bddDomain);
    end
    
    function W = setupGreensFunction(W)
        % Now for the Green's function.
        
        if W.theDomain.m > 1
            W.greensBeta = greensC0(W.beta, skpDomain(W.bddDomain));
        end
    end
end

end
