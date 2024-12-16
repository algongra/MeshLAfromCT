function [cux cuy cuz cp n indxrb] = genrefcent(dom0,domf,domr0,domrf,k0,kf,...
                                                ratio0,ratiof,dr,dir,varargin)
% @author A.Gonzalo
% 
% @brief function to build a mesh with:
%        coarse|stretching|refinement|stretching|coarse sections
% 
% @date  25-01-2015 by A.Gonzalo
%        Created
% @date  26-01-2015 by A.Gonzalo\n
%        Documented
% @date  09-04-2015 by A.Gonzalo\n
%        Added indxrb
% 
% @details
%
% @verbatim
% + Example of the mesh in y direction:
%
%  ------------------------------------------------------------ yf
%                      COARSE SECTION (dyc)
%  ------------------------------------------------------------
%            FINAL STRETCHING SECTION (dyr*(1+kyf)^i)
%  ------------------------------------------------------------ yrf
%                   FINAL REFINED SECTION (dyr)
%  ------------------------------------------------------------ yr0
%           INITIAL STRETCHING SECTION (dyr*(1+ky0)^i)
%  ------------------------------------------------------------
%                      COARSE SECTION (dyc)
%  ------------------------------------------------------------ y0
% @endverbatim
%
%  MANDATORY ARGUMENTS:
%  --------------------
%
%  - dom0: intial point of physical domain
%  - domf: final point of physical domain
%
%  - domr0: intial physical point of refined section
%  - domrf: final physical point of refined section
%
%  - k0: stretching factor of initial stretching section
%  - kf: stretching factor of final stretching section
%
%  - ratio0: refinement ratio of initial coarse section
%  - ratiof: refinement ratio of final coarse section
%
%  - dr: distance between points in refined section
%
%  - dir: direction you want to mesh:
%         x --> 1
%         y --> 2
%         z --> 3
%
%  OUTPUT ARGUMENTS:
%  -----------------
%
%  - cux: points of ux
%  - cuy: points of uy
%  - cuz: points of uz
%  - cp: points of p
%  - indxrb: beggining of refined zone index
%
% @code
% dir = 1;
% [xux2 xuy2 zuz2 xp2 nx irefbeg] = genrefcent(dom0,domf,domr0,domrf,k0,kf,...
%                                              ratio0,ratiof,dr,dir);
% dir = 2;
% [yux2 yuy2 yuz2 yp2 ny jrefbeg] = genrefcent(dom0,domf,domr0,domrf,k0,kf,...
%                                              ratio0,ratiof,dr,dir,...
%                                              'ghost',true);
% dir = 3;
% [zux2 zuy2 zuz2 zp2 nz krefbeg] = genrefcent(dom0,domf,domr0,domrf,k0,kf,...
%                                              ratio0,ratiof,dr,dir,...
%                                              'ghost',true,'per',perz);
% @endcode

per = 0;
ghost = false;
bcfix = false;
misc.assigndefaults(varargin{:});

% NUMBER OF SECTORS
nsector = 5; % number of sectors

% NUMBER OF POSITIONS (BEGGINING AND END)
npos = 2; % number of positions (beggining and end)

% POSITIONS OF THE MESH
iibeg = 1; % beggining
iiend = 2; % end

% COARSE RESOLUTION
dc = zeros(1,npos); % coarse resolution coordinate array

% STRETCHING RATIO
sr = zeros(1,npos); % stretching ratio coordinate array

% STRETCHING ZONE SIZE
ns = zeros(1,npos); % reference stretching zone size coordinate array

% PHYSICAL DOMAIN
pos = zeros(1,npos); % bounds coordinate array (physical domain)

% PHYSICAL DOMAIN OF SOLVER MESH
possp = zeros(1,npos); % bounds coordinate array (physical domain solver mesh)

% REFINED PHYSICAL DOMAIN
rpos = zeros(1,npos); % refined bounds coordinate array

% STRETCHING FACTOR 
k = zeros(1,npos); % stretching factor coordinate array

% REFINEMENT RATIO
ratio = zeros(1,npos); % refinement ratio coordinate array

% ARRAY USED TO SET POSITION OF THE MESH
iiarray1 = zeros(1,npos); % array used to set position of the mesh (iibeg,iiend)

% STRETCHING PHYSICAL DOMAIN
spos =zeros(1,npos); % stretching bounds

% COARSE ZONE SIZE
nc = zeros(1,npos); % reference coarse zone size


% Initialize the main arrays
% PHYSICAL DOMAIN
pos(iibeg) = dom0;
pos(iiend) = domf;

% REFINED PHYSICAL DOMAIN
rpos(iibeg) = domr0;
rpos(iiend) = domrf;

% STRETCHING FACTOR 
k(iibeg) = k0;
k(iiend) = kf;

% REFINEMENT RATIO
ratio(iibeg) = ratio0;
ratio(iiend) = ratiof;

% CALCULATE COARSE RESOLUTION
dc0 = dr*ratio0;
dcf = dr*ratiof;
dc(iibeg) = dc0;
dc(iiend) = dcf;

% ARRAY USED TO SET POSITION OF THE MESH
iiarray1 = [iibeg iiend];

% Calculate nr
% Calculate refined zone size
nr = (rpos(iiend)-rpos(iibeg))/dr;
nr = int32(nr);

for i2=1:npos
    iipos = iiarray1(i2);
    % Calculate sr
    % Calculate stretching ratio
    sr(iipos) = 1 + k(iipos);

    % Get ncs(iic,iipos) and leng
    % Get reference stretching zone size and length of stretching
    %  zone
    lstop = abs(rpos(iipos) - pos(iipos));
    dstop = dc(iipos);
    [nsaux leng] = fields.strnumpts(lstop,dstop,dr,sr(iipos));
    ns(iipos) = nsaux;

    % Calculate spos(iipos)
    % Calculate bounds of stretching sections (left and right)
    if iipos==iibeg
       % Calculate left bound of stretching section
       spos(iipos) = rpos(iipos) - leng;
    else
       % Calculate right bound of stretching section
       spos(iipos) = rpos(iipos) + leng;
    end

    if abs(spos(iipos)-pos(iipos)) >= dc(iipos)
       if iipos==iibeg
          % Get aproximate possp(iibeg) to calculate nc(iibeg)
          possp(iipos) = pos(iipos) + dc(iipos);

          % Calculate value of nc(iibeg)
          nc(iipos) = round((spos(iipos)-possp(iipos))/dc(iipos));

          % Get possp(iibeg)
          % Get initial coordinate point of sp mesh
          possp(iipos) = spos(iipos) - dc(iipos)*nc(iipos);

          % Get pos(iibeg)
          % Get new initial coordinate point of gp mesh
          pos(iipos) = possp(iipos) - dc(iipos);
       elseif iipos==iiend
          % Get aproximate possp(iiend) to calculate nc(iiend)
          possp(iipos) = pos(iipos) - dc(iipos);

          % Calculate value of nc(iiend)
          nc(iipos) = round((possp(iipos)-spos(iipos))/dc(iipos));

          % Get possp(iiend)
          % Get final coordinate point of sp mesh
          possp(iipos) = spos(iipos) + dc(iipos)*nc(iipos);

          % Get pos(iiend)
          % Get new final coordinate point of gp mesh
          pos(iipos) = possp(iipos) + dc(iipos);
       end
    else
       % Get the distance between gp mesh and sp mesh
       dc(iipos) = dr*sr(iipos).^ns(iipos);

       % Get possp(:)
       % Get initial and final coordinate point of sp mesh
       if iipos==iibeg
          % Get possp(iibeg)
          % Get initial coordinate point of sp mesh
          possp(iipos) = rpos(iipos) - leng + dc(iipos);

          if bcfix
             % Get pos(iibeg)
             % Get new initial coordinate point of gp mesh
             pos(iipos) = rpos(iipos) - leng;
          else
             % Get dc(iibeg)
             % Get the coarse size of the mesh at the beggining of the mesh
             dc(iipos) = dr*sr(iipos).^(ns(iipos)-1);

             % Get pos(iibeg)
             % Get new initial coordinate point of gp mesh
             pos(iipos) = possp(iipos) - dc(iipos);
          end
       elseif iipos==iiend
          % Get possp(iiend)
          % Get final coordinate point of sp mesh
          possp(iipos) = rpos(iipos) + leng - dc(iipos);

          if bcfix
             % Get pos(iiend)
             % Get new final coordinate point of gp mesh
             pos(iipos) = rpos(iipos) + leng;
          else
             % Get dc(iiend)
             % Get the coarse size of the mesh at the end of the mesh
             dc(iipos) = dr*sr(iipos).^(ns(iipos)-1);

             % Get pos(iiend)
             % Get new final coordinate point of gp mesh
             pos(iipos) = possp(iipos) + dc(iipos);
          end
       end
       ns(iipos) = ns(iipos) - 1;

       nc(iipos) = 0;
    end
end

% Calculate ic0b, ic0e, icfb, icfe,
%                 irb , ire ,
%           is0b, is0e, isfb, isfe
ic0b = 1;
ic0e = 1+nc(iibeg);
%is0b = 1+nc(iibeg)+1;
is0e = 1+nc(iibeg)+ns(iibeg);
indxrb  = 1+nc(iibeg)+ns(iibeg)+1;
ire  = 1+nc(iibeg)+ns(iibeg)+nr;
%isfb = 1+nc(iibeg)+ns(iibeg)+nr+1;
isfe = 1+nc(iibeg)+ns(iibeg)+nr+ns(iiend);
%icfb = 1+nc(iibeg)+ns(iibeg)+nr+ns(iiend)+1;
icfe = 1+nc(iibeg)+ns(iibeg)+nr+ns(iiend)+nc(iiend);

% Calculate n of body centered mesh
% Calculate solver point mesh size of body centered mesh
n = icfe + 1;

% build points of solver mesh [face centered mesh and body centered mesh]
u(ic0b) = possp(iibeg);
i = ic0b;
while ic0e > i
      i = i + 1;
      u(i) = u(i-1) + dc(iibeg);
      p(i) = (u(i) + u(i-1))/2;
end
j = ns(iibeg);
while is0e > i
      i = i + 1;
      u(i) = u(i-1) + dr*sr(iibeg)^j;
      p(i) = (u(i) + u(i-1))/2;
      j = j - 1;
end
while ire > i
      i = i + 1;
      u(i) = u(i-1) + dr;
      p(i) = (u(i) + u(i-1))/2;
end
j = 1; 
while isfe > i
      i = i + 1;
      u(i) = u(i-1) + dr*sr(iiend)^j;
      p(i) = (u(i) + u(i-1))/2;
      j = j + 1;
end
while icfe > i
      i = i + 1;
      u(i) = u(i-1) + dc(iiend);
      p(i) = (u(i) + u(i-1))/2;
end

% finish building points of solver mesh [body centered mesh]
p(ic0b) = p(ic0b+1) - dc(iibeg);
p(n) = p(n-1) + dc(iiend);

% build points of ghost mesh [face centered mesh]
u0 = zeros(1,icfe+2);
u0(ic0b+1:icfe+1) = u;
u0(ic0b) = u(ic0b) - dc(iibeg);
u0(icfe+2) = u(icfe) + dc(iiend);

% build points of ghost mesh [body centered mesh]
p0 = zeros(1,n+2);
p0(ic0b+1:n+1) = p;
if bcfix
   p0(ic0b) = p(ic0b) - dc(iibeg)*sr(iibeg);
   p0(n+2) = p(n) + dc(iiend)*sr(iiend);
else
   p0(ic0b) = p(ic0b) - dc(iibeg);
   p0(n+2) = p(n) + dc(iiend);
end

% add one point if periodic
if per == 1
   % recall n = icfe + 1
   u(n) = u(icfe) + dc(iiend);
   u0(n+2) = u0(icfe+2) + dc(iiend);
end
 
% check periodicity of the mesh
if per == 1 && dc(iibeg) ~= dc(iiend)
   funcnm = mfilename;
   disp(' ')
   disp(' ')
   fprintf('ERROR SETTING MESH PARAMETERS IN %s FUNCTION',funcnm)
   disp(' ')
   disp(' ')
   aggfrmt = ['A periodic mesh must have the same coarse spacing at the'...
              'beggining and at the end [dc(iibeg) = dc(iiend)]'];
   disp(aggfrmt)
   quit
end

% define output arrays
if dir == 1
   cux = u;
   cuy = p;
   cuz = p;
   cp = p;
   indxrb = [indxrb;indxrb+1;indxrb+1];
   
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
   indxrb = [indxrb+1;indxrb;indxrb+1];
   
   if ghost
      cux = p0;
      cuy = u0;
      cuz = p0;
      cp = p0;
   end
elseif dir == 3
   cux = p;
   cuy = p;
   cuz = u;
   cp = p;
   indxrb = [indxrb+1;indxrb+1;indxrb];
   
   if ghost
      cux = p0;
      cuy = p0;
      cuz = u0;
      cp = p0;
   end
end

return
end
