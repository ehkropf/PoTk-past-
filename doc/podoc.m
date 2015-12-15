function podoc(varargin)
%PODOC publishes the PoTk document set.

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

dlist = doclist;

mxdom = 'mxdom2html.xsl';
exargs = {};

if nargin 
    % Is first argument a doc file name? If so, process just that file.
    try
        vstr = validatestring(varargin{1}, dlist);
    catch
        vstr = [];
    end
    if ~isempty(vstr)
        match = strcmp(vstr, dlist(:,1));
        dlist = dlist(match, :);
        if nargin > 1
            varargin = varargin(2:end);
        else
            varargin = {};
        end
    end
    
    % Option stuff.
    if numel(varargin) > 0
        switch varargin{1}
            case 'openonly'
                if exist(mname2htname(dlist{1,1}), 'file')
                    openup(dlist)
                    return
                end
                
            case 'textonly'
                exargs = {'evalCode', false};
        end
    end
end

hadupdir = any(strcmp(strsplit(path, pathsep), '..'));
if ~hadupdir
    addpath('..')
end

if size(dlist, 1) > 0
    fprintf(['Building documentation. This may take a bit.\n' ...
        'Please be patient.\n'])
    tic
    try
        for k = 1:size(dlist, 1)
            srcfile = dlist{k,1};
            publish(srcfile, 'stylesheet', mxdom, exargs{:});
        end
    catch error
        toc
        cleanup(hadupdir)
        rethrow(error)
    end
    
    toc
    cleanup(hadupdir)
    openup(dlist)
end

end

%--------------------------------------------------------------------
function cleanup(hadupdir)

if ~hadupdir
    rmpath('..')
end

end

function openup(dlist)

% Assume doc to open is first in the list.
toopen = mname2htname(dlist{1,1});

try
    open(toopen)
catch err
    error('PoTk:RuntimeError', ...
        'Unable to open document file "%s".\n%s', ...
        toopen, err.message)
end

end

function fname = mname2htname(fname)

[~, fname, ~] = fileparts(fname);
fname = ['html/' fname '.html'];    

end
