%% Basic multiply connected flow examples.

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


%%
% Some setup for use below.

clear
aspectequal = @(ax) set(ax, 'dataaspectratio', [1, 1, 1]);
streamColor = [0, 0.447, 0.741];
vectorColor = [0.929, 0.694, 0.125];


%%
% Some bounded flow. Pick some circles and some points.

dv = [
    -0.25524+0.38986i
    -0.08042-0.18706i
    0.54196-0.30944i];
qv = [
    0.12362
    0.15385
    0.12607];
