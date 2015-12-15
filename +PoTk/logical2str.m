function str = logical2str(tfval)
%LOGICAL2STR converts logical to string.

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

% if numel(tfval) == 1
%     if tfval
%         str = 'true';
%     else
%         str = 'false';
%     end
%     return
% end

str = cell(size(tfval));
str(:) = {'false'};
str(tfval) = {'true'};
if numel(str) == 1
    str = str{:};
end
