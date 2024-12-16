function [xf varargout] = mesher(xy,dx,varargin)
% @author M.Moriche
%
% @brief
%
% @date 03-11-2014 M.Moriche \n
%       Created from airfoil_4d
%
% @details
%
% OUTPUT
%   - xf: coordiantes of the discretized geometry
%
% OPTIONAL OUTPUT
%   - ds: actual distance between consecutive points (archlength, which
%         means that if the discretization is poor, the distance between
%         points might be greater than this value.
%   - vol: volume associated to each point. (Lagrangian mesh for Immersed
%          Boundary method related)
%  
% MANDATORY ARGUMENTS
%   - xy: xy(i,j) coordinate of dimension i of point j
%   - dx: distance between consecutive points in discretized profile [double].
%
% OPTIONAL ARGUMENTS
%   - closed = true: Flat to generate closed curves.
%   - maxangle: Maximum angle permited when evolving along the chord of
%               the airfoil (inradians) [double].
%   - maxds: Maximum step between consecutive points. [double]
%
% To set variable arguments add as argument the name of the argument as string
% and the value of the argument. See examples below:
% 
%
% EXAMPLES:
% @code
% xf = mesher(x, y, dx, 'maxangle',0.001)
% 
% [xf ds] = mesher(xy, dx)
% [xf ds vol] = mesher(xy, dx)
% @endcode

% DEFAULTS
closed=true;
misc.assigndefaults(varargin{:});

n = size(xy,1);
p0=zeros(n,1);
perimeter = 0;
for i = 2:n
   perimeter = perimeter + norm(xy(i,:)-xy(i-1,:));
   p0(i,:) = perimeter;
end
perimeterref = perimeter;

if closed
   nreal = round(perimeterref/dx);
else
   nreal = round(perimeterref/dx)+1;
end

ds    = perimeter/nreal;
vol   = ds*dx;
xf=zeros(nreal,2);
perimeter = 0;
ix = 2;
xf(1,:) = xy(1,:);
dsleft = ds;
i=2;
while perimeter < perimeterref
    dsi = min(norm(xy(i,:)-xy(i-1,:)),perimeterref-perimeter);
    perimeternew = perimeter + dsi;

    %if perimeternew-perimeter >= dsleft
    if perimeternew-perimeter >= dsleft & ix <= nreal
       perimeter = perimeter + dsleft;
       x = interp1(p0, xy(:,1), perimeter);
       y = interp1(p0, xy(:,2), perimeter);
       xf(ix,:) = [x,y];
       ix = ix + 1;
       dsleft = ds;
    else
       perimeter = perimeternew;
       dsleft = dsleft - dsi;
       i = i +1;
    end
end
if ~closed
   xf(ix,:) = xy(end,:);
end

if nargout == 2
   varargout{1} = ds;
elseif nargout == 3
   varargout{1} = ds;
   varargout{2} = vol;
end

return 
end

