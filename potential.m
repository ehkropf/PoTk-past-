function pobj = potential(physicalDomain, varargin)
%POTENTIAL generates a complex potential object.

% Everett Kropf, 2015

if isa(physicalDomain, 'flowRegion')
    if isbounded(physicalDomain.islands)
        pobj = potentialInt(physicalDomain, varargin{:});
    else
        pobj = potentialExt(physicalDomain, varargin{:});
    end
else
    args = [physicalDomain, varargin(:)];
    pobj = builtin('potential', args{:});
end
