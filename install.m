function install(arg)
% PoTk install function.

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
