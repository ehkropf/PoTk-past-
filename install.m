function install(arg)
% PoTk install function.

% E. Kropf, 2014

if ~(exist('./potool', 'file') == 2 && ...
        exist('./potential', 'file') == 2 && ...
        exist('potential', 'class') == 8)
    error('PoTk:RuntimeError', ...
        'This utility must be run from the PoTk directory.')
end

if nargin && strcmp(arg, '-u')
    uninstall()
    return
end

fprintf('Checking to see if CMT is in your path ... ')
if exist('conformalmap', 'class') == 8 && exist('closedcurve', 'class') == 8
    fprintf('FOUND\n')
    ccpath = which('closedcurve.m');
    ccpath = fileparts(ccpath);
    fprintf('\tCMT found at path: %s\n', ccpath)
else
    fprintf('NO\n')

    cmtdir = [pwd '/CMT'];
    if exist(cmtdir, 'dir') == 7
        fprintf('Moving existing ./CMT to ./CMT.bak ... ')
        movefile('CMT', 'CMT.bak', 'f')
        fprintf('OK\n')
    end

    fprintf('Downloading CMT and unpacking ... ')
    % Do an actualy download here.
    unzip('https://github.com/ehkropf/conformalmapping/archive/master.zip')
    movefile('conformalmapping-master', 'CMT', 'f')
    fprintf('OK\n')
    
    fprintf('Adding ./CMT to your MATLAB path ... ')
    addpath(cmtdir);
    savepath
    fprintf('OK\n')
end

fprintf('PoTk install complete. Try ''podoc'' or ''potool'' now.\n')

end

function uninstall()

fprintf('Are we using a downladed CMT? ... ')
if exist('./CMT', 'dir') == 7
    fprintf('YES\n')
    
    fprintf('Removing ./CMT from your path ... ')
    cmtdir = [pwd '/CMT'];
    rmpath(cmtdir)
    savepath
    fprintf('OK\n')
else
    fprintf('NO\n')
end

fprintf('PoTk uninstall complete. This directory may be safely deleted.\n')

end
