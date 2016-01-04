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
        W.inputDomain = theDomain;
        
        W = constructPotential(W);
    end
end

methods(Access=protected)
    function w = calcPotential(W, z)
        
    end
    
        
    end
    
        
        end
    end
    
        
    end
    
        
        else
            zeta = mobius(1, 0, 0, 1);
        end
        
    end
    
        
        
        end
    end
    
%     function W = setupUniform(W)
%         % Precompute uniform flow components.
%         
%         if W.inputDomain.uniformStrength == 0 || ...
%                 W.inputDomain.m <= 1
%             return
%         end
% 
%         Chi = W.inputDomain.uniformAngle;
%         h = W.hCenterDiff;
%         
%         dbet = [...
%             W.beta + 0.5*h   % betaR
%             W.beta - 0.5*h   % betaL
%             W.beta + 0.5i*h  % betaU
%             W.beta - 0.5i*h  % betaD
%         ];
%         W.dbeta = dbet;
%         W.G0bxy = cell(2, 2);
%         
%         if mod(Chi, pi) > eps(pi)
%             % Sine term in use. Setup G0bx.
%             wb1 = skprime(dbet(1), W.wbeta);
%             wb1c = invParam(wb1);
%             wb2 = skprime(dbet(2), W.wbeta);
%             wb2c = invParam(wb2);
%             W.G0bxy(:,1) = { ...
%                 @(z) log(wb1(z)./wb1c(z)/abs(dbet(1)))/(2i*pi);
%                 @(z) log(wb2(z)./wb2c(z)/abs(dbet(2)))/(2i*pi)
%             };
%         end
%         if mod(Chi + pi/2, pi) > eps(pi)
%             % Cosine term in use. Setup G0by.
%             wb3 = skprime(dbet(3), W.wbeta);
%             wb3c = invParam(wb3);
%             wb4 = skprime(dbet(4), W.wbeta);
%             wb4c = invParam(wb4);
%             W.G0bxy(:,2) = { ...
%                 @(z) log(wb3(z)./wb3c(z)/abs(dbet(3)))/(2i*pi);
%                 @(z) log(wb4(z)./wb4c(z)/abs(dbet(4)))/(2i*pi)                
%             };
%         end
%     end
end

end
