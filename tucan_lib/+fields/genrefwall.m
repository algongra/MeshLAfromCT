function [cux cuy cuz cp dr dc indxrb] = genrefwall(dom0,domf,n,k,dir,varargin)

% @author A.Gonzalo
% 
% @brief function to build one direction of a 3D mesh with better
%        resolution near the bounds of the mesh
% 
% @date  20-01-2015 by A.Gonzalo\n
%        Created
% @date  20-01-2015 by A.Gonzalo\n
%        Documented
% @date  17-09-2015 by A.Gonzalo\n
%        Added indxrb
% 
% @details
%
% -----------------------------------------------------------------------------
% + YOU CAN NOT USE THIS FUNCTION TO CREATE THE EULERIAN MESH IF YOU ARE GOING
%   TO USE IMMERSED BOUNDARIES (EXPLANATION BELOW IN OUTPUT ARGUMENTS [indxrb])
% -----------------------------------------------------------------------------
% + This function uses a stretching factor to get better resolution near the
%   bounds of the mesh (e.g. In a channel, better resolution near the walls)
% + The mesh will be always symmetrical
%
% + MANDATORY ARGUMENTS:
%
%  - dom0: initial point of physical domain [ghost points of body mesh]
%  - domf: final point of physical domain [ghost points of body mesh]
%  - n:    number of point in the mesh [solver points of body mesh]
%  - k:    stretching factor
%  - dir:  direction (x=1, y=2, z=3)
%
% + OUTPUT ARGUMENTS:
%
%  - cux: x velocity points [face/body centered mesh]
%  - cuy: y velocity points [face/body centered mesh]
%  - cuz: z velocity points [face/body centered mesh]
%  - cp:  pressure points [body centered mesh]
%  - dr:  smallest distance between points (dxr, dyr, dzr)
%  - dc:  larger distance between points (dxr, dyr, dzr)
%  - indxrb: beggining of refined zone index.
%            -------------------------------------------------------------------
%            EXPLANATION:
%
%            IN THIS FUNCTION, SET THIS PARAMETER HAS NO MEANING.
%            THIS FUNCTION IS USEFUL TO SOLVE PROBLEMS WHERE YOU NEED MANY
%            POINTS NEAR THE WALLS OF THE EULERIAN MESH (CHANNELS), THEREFORE IT
%            CREATES TWO REFINED ZONES (ONE ON EACH WALL).
%            TUCAN3REF ONLY WORKS WITH IMMERSED BOUNDARIES WHEN THERE IS ONLY
%            ONE REFINED ZONE.
%            I'M ADDING THE PARAMETER BECAUSE IS A MANDATORY INPUT IN TUCAN3REF
%            MESHES.
%            -------------------------------------------------------------------
%
% @code
% dir=1;
% [xux xuy xuz xp dxr dxc irefbeg] = fields.genrefwall(dom0,domf,n,k,dir,'per',1)
%
% dir=3;
% [zux0 zuy0 zuz0 dzr dzc krefbeg] = fields.genrefwall(dom0,domf,n,k,dir,...
%                                                      'ghost',true,'per',1);
% @endcode

per = 0;
ghost = false;
misc.assigndefaults(varargin{:});

% warning
if mod(n,2) ~= 0
   fprintf('\nn is odd. The mesh must have a number of points (n) even\n\n')
   quit
end

% number of points of face centered mesh [solver points]
m = n-1;
N = n/2;
M = N-1;
% number of points of face centered mesh [ghost points]
mg = m+2;
Mg = M+2;

% get stretching factor
sf = 1 + k;

% build face centered mesh [ghost points]
u0(1) = 0; u0(2) = 1;
for i=3:Mg;
    du = (u0(i-1)-u0(i-2))*sf;
    u0(i) = u0(i-1) + du;
end

% reescale the mesh
reescale_f = ((domf-dom0)/2)/u0(Mg);
u0(1:Mg) = u0(1:Mg)*reescale_f;

% build symmetric points
j = 1;
for i=(Mg+1):mg;
    u0(i) = u0(i-1) + u0(i-j)-u0(i-j-1);
    j = j + 2;
end

% move the mesh to the initial point
u0(:) = u0(:) + dom0;

% build face centered mesh [solver points]
u = zeros(1,m);
u = u0(2:end-1);

% Get distance between points in the wall [needed for periodic and ghost meshes]
dr = u(2) - u(1);
% Get distance in the center of the channel
dc = max(diff(u));

% rebuild face centered mesh [ghost points]
u0 = zeros(1,mg);
u0(2:end-1) = u;
u0(1) = u(1) - dr;
u0(end) = u(end) + dr;

% build body centered mesh [solver points]
p = zeros(1,n);
for i=1:n;
    p(i) = (u0(i+1)+u0(i))/2;
end

% build body centered mesh [ghost points]
p0 = zeros(1,n+2);
p0(2:end-1) = p;
p0(1) = p(1) - dr;
p0(end) = p(end) + dr;

% add one point if periodic
if per == 1
   m = m + 1;
   u(m) = u(m-1) + dr;
   u0(m+2) = u0(m+1) + dr;
end

indxrb = zeros(3,1);

% define output arrays
if dir == 1
   cux = u;
   cuy = p;
   cuz = p;
   cp = p;
   
   if ghost
      cux = u0;
      cuy = p0;
      cuz = p0;
      cp = p0;
   end
elseif dir == 2
   cux = p;
   cuy = u;
   cuz = p;
   cp = p;
   
   if ghost
      cux = p0;
      cuy = u0;
      cuz = p0;
      cp = p0;
   end
else
   cux = p;
   cuy = p;
   cuz = u;
   cp = p;
   
   if ghost
      cux = p0;
      cuy = p0;
      cuz = u0;
      cp = p0;
   end
end

return
end
