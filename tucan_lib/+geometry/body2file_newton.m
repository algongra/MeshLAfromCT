function body2file(xyz,vol,xyzc,ib,ie,spi,rhoratio,bodyvol,fnm,varargin)
% @author M.Moriche
%
% @brief Function to write coordinates of Lagrangian mesh defined by
%        the matrix xy and associated volume vol in file fnm.
%        With different body handling
%
% @date 02-10-2014 by M.Moriche \n
%       Created from grid2file
% @date 21-05-2015 by M.Moriche \n
%       density ratios and particle volume included
%
% @details
%
% NO OUTPUT
%
% MANDATORY ARGUMENTS
%  - xyz(n,ndim): matrix with coordinates of Lagrangian grid points [double]
%  - vol(n): associated volume of each Lagrangian marker. [double]
%  - xyzc(nbody,ndim): control point for each body [double]
%  - ib(nbody): begin index of each body [integer]
%  - ie(nbody): ending index of each body [integer]
%  - spi(ntensor, nbody): specific inertia of each body
%  - rhoratio(nbody): rho_p/rho_f of each body
%  - bodyvol(nbody): volume of each body
%  - fnm: file name to write. [string]
%
% OPTIONAL ARGUMENTS
%  - prec: precision of real numbers. [integer] = 8
%  - path: path to save the file [string] = '.'
%
%  EXAMPLES:
%
%  @verbatim
%  body2file(xy,vol,xyc,ib,ie,spi,rhoratio, bodyvol, fnm)
%  body2file(xyz,vol,xyc,ib,ie,spi,fnm,rhoratio, bodyvol, varargin)
%  body2file(xy,vol,xyc,ib,ie,spi,fnm,rhoratio, bodyvol, 'path','/data/geometry')
%  @endverbatim
%

% defaults
prec = 8;          % real numbers precision
path = '.';

misc.assigndefaults(varargin{:});

[nreal ndim] = size(xyz);
[nbody ndim] = size(xyzc);

ntensor = 0.5*ndim*(ndim+1);

fullfnm = fullfile(path, fnm);
%
fid=fopen(fullfnm,'w','l');
% precision
fwrite(fid,4,   'int');
fwrite(fid,prec,'int');
fwrite(fid,4,   'int');
% dimensions
fwrite(fid,4,   'int');
fwrite(fid,ndim,'int');
fwrite(fid,4,   'int');
% lagrangian size
fwrite(fid,4, 'int');
fwrite(fid,nreal, 'int');
fwrite(fid,4, 'int');
% lagrangian mesh
for i=1:ndim; 
   fwrite(fid,8*nreal, 'int');
   fwrite(fid,xyz(:,i),   'double');
   fwrite(fid,8*nreal, 'int');
end
% marker volume
fwrite(fid,8*nreal,   'int');
fwrite(fid,vol, 'double');
fwrite(fid,8*nreal,   'int');
% number of bodies
nbody
fwrite(fid,4, 'int');
fwrite(fid,nbody, 'int');
fwrite(fid,4, 'int');
% control point
for i=1:ndim 
   fwrite(fid,8*nbody,     'int');
   fwrite(fid,xyzc(:,i),'double');
   fwrite(fid,8*nbody,     'int');
end 
% begin index
fwrite(fid,4*nbody,'int');
fwrite(fid,ib,'int');
fwrite(fid,4*nbody,'int');
% end index
fwrite(fid,4*nbody,'int');
fwrite(fid,ie,'int');
fwrite(fid,4*nbody,'int');
% specific inertia (I/M)
fwrite(fid,8*nbody*ntensor,'int');
fwrite(fid,spi,'double');
fwrite(fid,8*nbody*ntensor,'int');
% density ratio
fwrite(fid,8*nbody        ,'int');
fwrite(fid,rhoratio,    'double');
fwrite(fid,8*nbody        ,'int');
% body volume
fwrite(fid,8*nbody        ,'int');
fwrite(fid,bodyvol,     'double');
fwrite(fid,8*nbody        ,'int');

fclose(fid);

return 
end

