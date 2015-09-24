function str = logical2str(tfval)
%LOGICAL2STR converts logical to string.

% E. Kropf, 2014

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
