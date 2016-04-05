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
    zetaFun                 % (r/o) Mobius transform unbounded -> bounded.
    beta = inf              % (r/o) Image of infinity in bounded domain.
    greensBeta              % (r/o) Green's function wrt C0 at beta.
end

properties(Access=protected)
    defaultPlotScale = 1.5
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
    
    function Wb = potentialBdd(W)
        %Convert external to bounded potential.
        
        Wb = W.bddPotential;
    end
end

methods(Access=protected)
    %%%%% Calculation.
    function w = calcPotential(W, z)
        % Combine with internal calculation.
        
        W = waitbarInitialize(W, 'Computing potential values');
        zeta = W.zetaFun(z);
        w = calcBounded(W, zeta);
        w = w - calcExternalPart(W, zeta);
        waitbarRelease(W);
    end
    
    function w = calcBounded(W, zeta)
        % Bounded potential part.
        
        if W.useWaitBar
            subbar = PoG.subbar(W.waitBar, 3/4);
            w = calcPotential(W.bddPotential, zeta, subbar);
        else
            W.bddPotential.useWaitBar = W.useWaitBar;
            w = W.bddPotential(zeta);
        end
    end
    
    function w = calcExternalPart(W, zeta)
        % Calculate extra part for external flow.
        
        waitbarUpdate(W, 3/4, 'Externalizing flow')
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
        
        W = waitbarInitialize(W, 'Configuring potential functions');
        W = setupBoundedRegion(W);
        W = setupBoundedPotential(W);
        W = setupGreensFunction(W);
        W = waitbarRelease(W);
    end
    
    function W = setupBoundedRegion(W)
        % Mobius transform of exterior domain.
        
        waitbarUpdate(W, 0, 'Equivalent bounded domain')
        if W.theDomain.m > 0
            [Db, zeta] = regionBdd(W.theDomain);
            
            W.beta = pole(inv(zeta));
            W.zetaFun = zeta;
            W.bddDomain = Db;
        else
            % No equivalent bounded domain.
            W.zetaFun = mobius(1, 0, 0, 1);
        end
    end
    
    function W = setupBoundedPotential(W)
        % Get the bounded potential.
        
        waitbarUpdate(W, 1/5, 'Prime domain potential')
        nv = numel(W.varPropList);
        args = cell(1, 2*nv);
        for i = 1:nv
            var = W.varPropList{i};
            args(2*i-[1,0]) = {var, W.(var)};
        end
        if W.useWaitBar
            args = [{PoG.subbar(W.waitBar, 3/5)}, args];
        end
        W.bddPotential = potentialBdd(W.bddDomain, args{:});
    end
    
    function W = setupGreensFunction(W)
        % Now for the Green's function.
        
        waitbarUpdate(W, 4/5, 'Singularity at infinity')
        if W.theDomain.m > 1
            W.greensBeta = greensC0(W.beta, skpDomain(W.bddDomain));
        end
    end
end

end
