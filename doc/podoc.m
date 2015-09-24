function podoc(varargin)
%PODOC publishes the PoTk document set.

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
