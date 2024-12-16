function body2fileS(lg,fnm,varargin)
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
% @date 01-06-2016 by A.Gonzalo \n
%       normal vectors of each lagrangian point included as optional argument
% @date 23-11-2020 by E.Duran \n
%       input arguments in struct format
% @date 09-02-2022 by A.Gonzalo \n
%       setting back varargin to define the path and prec of the file
%
% @details
%
% NO OUTPUT
%
% MANDATORY ARGUMENTS
%  + lg: variable containing the following fields [struct]:
%    - xyz(n,ndim): matrix with coordinates of Lagrangian grid points [double]
%    - vol(n): associated volume of each Lagrangian marker. [double]
%    - xyzc(nbody,ndim): control point for each body [double]
%    - normvec(n,ndim): normal vector of each Lagrangian point [double]
%    - ib(nbody): begin index of each body [integer]
%    - ie(nbody): ending index of each body [integer]
%    - spi(ntensor, nbody): specific inertia of each body
%  + fnm: file name to write. [string]
%
% OPTIONAL ARGUMENTS
%  + prec: precision of real numbers. [integer] = 8
%  + path: path to save the file [string] = '.'
%
%  EXAMPLES:
%
%  @verbatim
%  body2fileS(lg,fnm)
%  body2fileS(lg,fnm,varargin)
%  body2fileS(lg,spi,fnm,'path','/data/geometry')
%  @endverbatim


% defaults
prec = 8; % real numbers precision
path = '.';
misc.assigndefaults(varargin{:});

% Rename the information from lg struct
xyz = lg.xyz;
vol = lg.vol;
xyzc = lg.xyzc;
normvec = lg.normvec;
ib = lg.ib;
ie = lg.ie;
spi = lg.spi;
rhoratio = lg.rhoratio;
bodyvol = lg.bodyvol;
prec = lg.prec;

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
% normal vectors
for i=1:ndim
    fwrite(fid,8*nreal,     'int');
    fwrite(fid,normvec(:,i),'double');
    fwrite(fid,8*nreal,     'int');
end
% number of bodies
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

