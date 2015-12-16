classdef baseRegion
%baseRegion is PoTk base circle region class.

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
    centers                     % Circle centers column vector
    radii                       % Circle radii column vector
    circulation                 % Vector of scalars, boundary circulation
    singularities               % Point singularity location vector
    singStrength                % Vector of scalars for singularity strength
    uniformStrength = 0         % Uniform background field strength
    uniformAngle = 0            % Background field angle
end

properties(Dependent, Hidden)
    m                           % Connectivity
end

methods
    function R = baseRegion(dv, qv, varargin)
        if ~nargin
            return
        end
        
        switch class(dv)
            case 'skpDomain'
                varargin = [qv, varargin];
                [dv, qv] = domainData(dv);
                
            case 'circleRegion'
                varargin = [qv, varargin];
                qv = dv.radii;
                dv = dv.centers;
                
            case 'double'
                dv = dv(:);
                qv = qv(:);
                if ~isequal(size(dv), size(qv))
                    error(PoTk.ErrorTypeString.InvalidArgument, ...
                        'Center and radius vectors must be the same size.')
                end
                
            otherwise
                error(PoTk.ErrorTypeString.InvalidArgument, ...
                    ['Expected first argument to be a vector of '...
                    'centers, an skpDomain object, or a circleRegion '...
                    'object.'])
        end
        
        R.centers = dv;
        R.radii = qv;
        
        nvararg = numel(varargin);
        if nvararg == 0
            R.circulation = zeros(size(dv));
        else
            props = properties('baseRegion');
            props = props(3:end);
            if ischar(varargin{1})
                for i = 1:2:nvararg-1
                    try
                        vstr = validatestring(varargin{i}, props);
                        R.(vstr) = varargin{i+1};
                    catch err
                        switch err.identifier
                            case 'MATLAB:unrecognizedStringChoice'
                                errmsg = sprintf(...
                                    'Parameter ''%s'' unrecognized.', ...
                                    varargin{i});
                            case 'MATLAB:badsubsript'
                                errmsg = sprintf(...
                                    'Value expected for parameter ''%s''.', ...
                                    varargin{i});
                            otherwise
                                rethrow(err)
                        end
                        error(PoTk.ErrorTypeString.InvalidArgument, ...
                            '%s', errmsg)
                    end
                end
            else
                for i = 1:nvararg
                    try
                        R.(props{i}) = varargin{i};
                    catch err
                        error(PoTk.ErrorTypeString.InvalidArgument, ...
                            '%s \nArguments not recognized. See help.', ...
                            err.message)
                    end
                end
            end
        end
    end
    
    %---------------------------------------------------------------------
    function m = get.m(R)
        m = mGetter(R);
    end
    
    function R = set.circulation(R, vec)
        if any(~PoTk.isRealAndFinite(vec))
            error(PoTk.ErrorTypeString.InvalidValue, ...
                'Circulation values must be real and finite.')
        end
        R.circulation = vec(:);
    end
    
    function R = set.singularities(R, vec)
        if any(~PoTk.isFiniteNumber(vec))
            error(PoTk.ErrorTypeString.InvalidValue, ...
                'Singularity locations must be finite.')
        end
        R.singularities = vec(:);
    end
    
    function R = set.singStrength(R, vec)
        if any(~PoTk.isRealAndFinite(vec))
            error(PoTk.ErrorTypeString.InvalidValue, ...
                'Singularity strength values must be real and finite.')
        end
        R.singStrength = vec(:);
    end
end

methods(Access=protected,Abstract)    
    m = mGetter(R)      % Return valid value for connectivity m.
end

methods(Access=protected)
    function sanityCheck(R)
        %Checks that region is in a valid state. Throws error if not.
        
        if ~isequal(numel(R.circulation), numel(R.centers))
            error(PoTk.ErrorTypeString.InvalidValue, ...
                ['Circulation vector must have the same number of\n'...
                'elements as there are inner boundaries.'])
        end
    end
end

end
