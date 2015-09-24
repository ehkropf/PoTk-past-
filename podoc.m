function podoc(docname)
%PODOC opens the PoTk documentation.
%
% podoc userguide
%   Opens the user guide.

% E. Kropf, 2014

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
