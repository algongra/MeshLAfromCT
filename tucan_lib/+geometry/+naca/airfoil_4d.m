function [xf varargout] = airfoil_4d(digits, c, dx, varargin)
% @author M.Moriche
%
% @brief Discretize NACA airfoils of 4 digits type for Immersed Boundary
%        Lagrangian grids.
%
% @date Far away in time by M.Moriche 
%       Created 
% @date November 2013 by M.Moriche & E.Hurtado
%       Review and bugs fixing
% @date November 2013 by E.Hurtado
%       Organized and documented.
% @date 3 December 2013 by M.Moriche
%       Adaptative step along chord. Move script to function.
% @date 13-12-2013 M.Moriche \n
%       Documented and organized with the new package files
% @date 30-12-2014 M.Moriche \n
%       Bug in first point (typo zerothicknes instead of zerothickness) \n
%       Path to checkstep was wrong
% @date 21-09-2015 M.Moriche \n
%       Bug in first last point computation for some dx,c realtion in zero thicnkess
%
% @details
%
% OUTPUT
%   - xf: coordiantes of the discretized airfoil.
%
% OPTIONAL OUTPUT
%   - ds: actual distance between consecutive points (archlength, which
%         means that if the discretization is poor, the distance between
%         points might be greater than this value.
%   - vol: volume associated to each point. (Lagrangian mesh for Immersed
%          Boundary method related)
%  
% MANDATORY ARGUMENTS
%   - digits: digits of the NACA proile [string]
%   - c: chord of the profile [double]
%   - dx: distance between consecutive points in discretized profile [double].
%
% OPTIONAL ARGUMENTS
%   - alfadeg: angle of attack (in degrees) [double]
%   - zerothickness: Flag to mofidy last coefficient of NACA law in order
%                    to have zero thickness at the trailing edge [boolean]
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
% xf = airfoil_4d(digits, c, dx,'zerothickness',true)
% xf = airfoil_4d(digits, c, dx,'alphadeg',20,'zerothickness',true)
% xf = airfoil_4d(digits, c, dx,'maxangle',0.001)
% 
% [xf ds] = airfoil_4d(digits, c, dx)
% [xf ds vol] = airfoil_4d(digits, c, dx)
% @endcode

% DEFAULTS
alfadeg = 0;
zerothickness=true;
maxds = -1;
maxangle = 0.005;

misc.assigndefaults(varargin{:});

if maxds == -1
   maxds = dx/2;
end

% code___________________________________________________
%
alfa = -alfadeg*pi/180;   % minus for naca function reference
%
n0=100000;  % max value for size of x0.
x0=zeros(n0,2);  % Unless otherwise stated, we initialize our points at 0

% initial values
[xl xu yl yu]=geometry.naca.naca_4d(digits,c,c,'alpha',alfa,'zerothickness',zerothickness);
xnow = [xu,yu];
i=1;
x0(i,:) = [xu, yu];
i = i+1;

dxi = 10*dx; % initialized with the step like mesh width
cdone = 0;
while cdone<c
    flag = false; 
    n2 = 0;
    dxi = dxi*2;
    while ~flag
       dxi = min(dxi/2,c-cdone);

       x1=cdone+dxi;
       [xl xu yl yu]=geometry.naca.naca_4d(digits,c-x1,c,'alpha',alfa,'zerothickness',zerothickness);
       xnext_1 = [xu, yu];

       x2=cdone+dxi/2;
       [xl xu yl yu]=geometry.naca.naca_4d(digits,c-x2,c,'alpha',alfa,'zerothickness',zerothickness);
       xnext_2 = [xu, yu];

       flag = geometry.naca.checkstep(xnow,xnext_1,xnext_2,maxds,maxangle);

       n2 = n2 + 1;
    end

    xnow = xnext_1;
    x0(i,:) = xnow;
    i = i+1;
    cdone = cdone+dxi;

    if n2 == 1
      dxi = dxi*4;
    end
end
% backwards
cleft=c;
while cleft>0
    flag = false; 
    n2 = 0;
    dxi = dxi*2;
    while ~flag
       dxi = min(dxi/2,cleft);

       x1=cleft-dxi;
       [xl xu yl yu]=geometry.naca.naca_4d(digits,c-x1,c,'alpha',alfa,'zerothickness',zerothickness);
       xnext_1 = [xl, yl];

       x2=cleft-dxi/2;
       [xl xu yl yu]=geometry.naca.naca_4d(digits,c-x2,c,'alpha',alfa,'zerothickness',zerothickness);
       xnext_2 = [xl, yl];

       flag = geometry.naca.checkstep(xnow,xnext_1,xnext_2,maxds,maxangle);

       n2 = n2 + 1;
    end

    xnow = xnext_1;
    x0(i,:) = [xl, yl];
    i = i+1;

    cleft = cleft-dxi;
    if n2 == 1
      dxi = dxi*4;
    end
end
[xl xu yl yu]=geometry.naca.naca_4d(digits,c,c,'alpha',alfa,'zerothickness',zerothickness);
xnow = [xl, yl];
x0(i,:) = xnow;

n = i;
x0 = x0(1:n,:);
p0=zeros(n,1);
perimeter = 0;
for i = 2:n
   perimeter = perimeter + norm(x0(i,:)-x0(i-1,:));
   p0(i,:) = perimeter;
end
perimeterref = perimeter;
%
if zerothickness
   nreal = round(perimeter/dx);
   %nreal = floor(perimeter/dx);
else
   nreal = round(perimeter/dx)+1;
   %nreal = floor(perimeter/dx)+1;
end
ds    = perimeter/nreal;
vol   = ds*dx;
xf=zeros(nreal,2);
perimeter = 0;
ix = 2;
xf(1,:) = x0(1,:);
dsleft = ds;
i=2;
while perimeter < perimeterref
    dsi = min(norm(x0(i,:)-x0(i-1,:)),perimeterref-perimeter);
    perimeternew = perimeter + dsi;

    %if perimeternew-perimeter >= dsleft
    if perimeternew-perimeter >= dsleft & ix <= nreal
       perimeter = perimeter + dsleft;
       x = interp1(p0, x0(:,1), perimeter);
       y = interp1(p0, x0(:,2), perimeter);
       xf(ix,:) = [x,y];
       ix = ix + 1;
       dsleft = ds;
    else
       perimeter = perimeternew;
       dsleft = dsleft - dsi;
       i = i +1;
    end
end

if ~zerothickness
   xf(ix,:) = x0(end,:);
end

if nargout == 2
   varargout{1} = ds;
elseif nargout == 3
   varargout{1} = ds;
   varargout{2} = vol;
end

return 
end

