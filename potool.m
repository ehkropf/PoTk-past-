function pog = potool(adomain)
%potool launches the PoTk interactive tool.
%
% See also the userguide in podoc.

% E. Kropf, 2014

if nargin
    pog = PoG.mainWin('run', adomain);
else
    pog = PoG.mainWin('run');
end

if ~nargout
    clear pog
end
