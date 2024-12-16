function sphere_txt2dat(h,r,fnm,varargin)
% @author A.Gonzalo
%
% @brief Function to get the coordinates points of a sphere from the .txt file
%        generated with points_on_sphere_2.cpp program, and to write coordinates
%        of Lagrangian sphere defined by the matrix xyz and associated volume
%        vol in file fnm.dat.
%
% @date 02-10-2014 by A.Gonzalo \n
% @date 01-07-2016 by M.Moriche \n
%       Update call to body2file with rhoratio, bodyvol
%
% @details
% MANDATORY ARGUMENTS
%  - h: refined elements size (dxr, dyr or dzr). [double]
%  - r: radius of the sphere. [double]
%  - fnm: file name (without .txt). [string]
%
% OPTIONAL ARGUMENTS
%  - centpos(ndim):  position of the sphere center. [double]
%  - prec: precision of real numbers. [integer] = 8
%  - path: path where .txt file is located [string] = '.'
%  - ifplot: plot sphere points. [false by default]
%
% MANDATORY ARGUMENTS TO CALL body2file.m
%  - xyz(n): matrix with coordinates of Lagrangian grid points. [double]
%  - vol(n): associated volume of each Lagrangian marker. [double]
%  - xyzc(nbody=1): control point for the body. [double]
%  - ib(nbody=1): begin index of each body. [integer]
%  - ie(nbody=1): ending index of each body. [integer]
%  - spi(ntensor,nbody=1): specific inertia of each body
%  - prec: precision of real numbers. [integer] = 8
%  - fnm: file name to write. [string]
%
% EXAMPLES:
%
% @code
% sphere_txt2dat(h,fnm)
% sphere_txt2dat(h,fnm,'ifplot',true)
% sphere_txt2dat(h,fnm,'path','/data/alex/IBcode/geometry')
% @endcode
%

% defaults
centpos = [0 0 0];
prec = 8;          % real numbers precision
path = '.';
ifplot = false;
misc.assigndefaults(varargin{:});


dottxt = '.txt';
dotdat = '.dat';

% READ data from file obtained by running the code points_on_sphere
fullfnm = strcat(fnm,dottxt);
fullfnm = fullfile(path,fullfnm);
[x, y, z] = textread(fullfnm,'%f %f %f');

% PARAMETERS NEEDED TO CALCULATE THE MANDATORY ARGUMENTS TO CALL body2file
% number of bodies
nbody = 1;
% number of dimensions
ndim  = 3;
% number of elements
nreal = length(x);

% CREATE MANDATORY ARGUMENTS REQUIRED TO CALL body2file
% lagrangian mesh
xyz = zeros(nreal,ndim);
xyz(:,1) = x*r + centpos(1);
xyz(:,2) = y*r + centpos(2);
xyz(:,3) = z*r + centpos(3);
% control point
xyzc = zeros(nbody,ndim);
% begin index of the body
ib(nbody) = 1;
% ending index of the body
ie(nbody) = nreal;
% marker volume
vol(1:nreal) = pi*h*(12*r^2+h^2)/(3*nreal);
% specific inertia (I/M)
ntensor = 0.5*ndim*(ndim+1);
spi = zeros(ntensor,nbody);

% clear useless arrays
clear x; clear y; clear z;

if ifplot
   % PLOT of the sphere and
   figure();
   plot3(xyz(:,1),xyz(:,2),xyz(:,3),'.'); axis equal; view(30,30);
end
% PRINT the dv
fprintf('\nvolume = %7.4e\n',vol(1));

% CALL body2file
fullfnm = strcat(fnm,dotdat)
%geometry.body2file(xyz,vol,xyzc,ib,ie,spi,fullfnm,'path',path);
geometry.body2file(xyz,vol,xyzc,ib,ie,spi,1,1,fullfnm,'path',path);

return 
end
