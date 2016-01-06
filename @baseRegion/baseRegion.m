classdef baseRegion
%baseRegion is PoTk base circle region class.
%
% This abstract class is not meant to be instantiated directly.
%
% Constructor:
% R = baseRegion(dv, qv)
%   Constructs basic region defined by circle center vector dv and radius
%   vector qv.
% R = baseRegion(dv, qv, circulation, singularities, singStrength, ...
%                uniformStrength, uniformAngle, dipole)
%   Constructs region by setting given parameters. The circulation vector
%   is bettter defined in the subclass. The singularities and singStrength
%   vectors must be the same size. The uniform field strength and angle are
%   real, finite values. The dipole is the location of the uniform field
%   dipole. Any parameter value may be empty per MATLAB function
%   convention.
% R = baseRegion(dv, qv, 'parameterName', parameterValue, ...)
%   This is the name/value pair version of the constructor.

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
    dipole = complex(inf)       % Uniform field dipole location
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
                        if isempty(varargin{i+1})
                            % Ignore empty argument.
                            continue
                        end
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
                        if isempty(varargin{i})
                            % Ignore empty argument.
                            continue
                        end
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
    
    function numberFeatures(R)
        %Number flow region features. Assumes plot has already been
        %called.
        
        numberBoundaries(R)
        d = axis;
        d = 0.01*max(diff(d(1:2)), diff(d(3:4)));
        for k = 1:numel(R.singularities)
            v = R.singularities(k);
            text(real(v) + d, imag(v) + d, num2str(k), 'color', 'm')
        end
    end
    
    function plot(R)
        %Draw the flow region.
        %
        % plot(R) where R is a baseRegion object.
        
        cah = newplot;
        washold = ishold;
        
        hold on
        Rc = circleRegion(R);
        fillHoles(R)
        plot(Rc)
        
        z = R.singularities;
        plot(real(z), imag(z), 'k.', 'markersize', 18)
        
        if ~washold
            cmtplot.whitefigure(cah)
            axis(plotbox(R))
            aspectequal
            axis off
            hold off
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

methods(Abstract)
    C = circleRegion(R)         % Convert to CMT circle region.
    numberBoundaries(R)         % Helper function for numberFeatures().
    xylim = plotbox(R, scale)   % Axis limits for plotting.
end

methods(Access=protected)
    function sanityCheck(R)
        % Call the sanity checks.
        
        try
            baseSanityCheck(R)
            subSanityCheck(R)
        catch err
            import PoTk.ErrorTypeString.*
            switch err.identifier
                case {InvalidValue, UndefinedState, InvalidArgument}
                    error(err.identifier, err.message)
                    
                otherwise
                    rethrow(err)
            end
        end
    end
    
    function baseSanityCheck(R)
        %Checks that region is in a valid state. Throws error if not.
        
        if numel(R.circulation) ~= numel(R.centers)
            error(PoTk.ErrorTypeString.InvalidValue, ...
                ['Circulation vector must have the same number of\n'...
                'elements as there are inner boundaries.'])
        end
        
        if numel(R.singularities) ~= numel(R.singStrength)
            error(PoTk.ErrorTypeString.InvalidValue, ...
                ['Singularity strength vector must have the same number\n'...
                'elements as the singularity location vector.'])
        end
        
        if isempty(R.uniformStrength) ...
                || ~PoTk.isRealAndFinite(R.uniformStrength)
            error(PoTk.ErrorTypeString.InvalidValue, ...
                'Uniform field strength must be a real, finite number.')
        end
    end
    
    function subSanityCheck(~)
        %Override this method in subclass for custom sanity checks.
    end
    
    function fillHoles(R)
        %Fill inner holes with fill color.
        
        dv = R.centers;
        qv = R.radii;
        for j = 1:numel(dv)
            fill(circle(dv(j), qv(j)))
        end
    end
end

methods(Access=protected,Abstract)
    m = mGetter(R)          % Return valid value for connectivity m.
end

end
