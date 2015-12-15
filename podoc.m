function podoc(docname)
%PODOC opens the PoTk documentation.
%
% podoc userguide
%   Opens the user guide.

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

if nargin
    docargs = {docname, 'openonly'};
else
    docargs = {'openonly'};
end

currDir = pwd;
potkDir = fileparts(which(mfilename));
if ~strcmp(potkDir, fileparts(which('flowRegion')))
    error('PoTk:runtimeError', 'Unable to locate PoTk document directory.')
end
cd([potkDir, '/doc'])

try
    podoc(docargs{:})
catch err
    cd(currDir)
    rethrow(err)
end
cd(currDir)
