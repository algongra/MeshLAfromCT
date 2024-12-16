function fou2fileS(fou,fnm,varargin)
% @author M.Garcia-Villalba
%
% @brief Function to write Fourier modes to calculate Lagrangian mesh points
%        coordinates and velocities at different time instants in a file 
%        readable for tucanH.
%
% @date ??-??-2018 by M.Garcia-Villalba \n
%       Created from body2file.m
% @date 23-11-2020 by E.Duran \n
%       input arguments in struct format
% @date 09-02-2022 by A.Gonzalo \n
%       added function documentation
%       setting back varargin to define the path and prec of the file
%
% MANDATORY ARGUMENTS
% hay que corregir estos tamaños de xyz, vecprod y separarlos en real/imaginario
%  + fou: variable containing the following fields [struct]:
%    - nreal: number of Lagrangian grid points (i.e., size(xyz_*,1)) [double]
%    - ndim: dimension of the agrangian grid (i.e., size(xyz_*,2)) [double]
%    - xyz_r(n,ndim,nm1): matrix with real part of Lagrangian grid points
%                         Fourier modes  [double]
%    - xyz_i(n,ndim,nm1): matrix with imaginary part of Lagrangian grid points
%                         Fourier modes [double]
%    - vecprod_r(n,ndim,nm2): matrix with real part of vector product
%                             Fourier modes [double]
%    - vecprod_i(n,ndim,nm2): matrix with imaginary part of vector product
%                             Fourier modes [double]
%    - nm1 : number of modes to be written for xyz
%    - nm2 : number of modes to be written for vecprod
%    - fnm : filename
%
% OPTIONAL ARGUMENTS
%  + prec: precision of real numbers. [integer] = 8
%  + path: path to save the file [string] = '.'
%
%  EXAMPLES:
%
%  @verbatim
%  fou2fileS(fou,fnm)
%  fou2fileS(fou,fnm,varargin)
%  fou2fileS(fou,spi,fnm,'path','/data/geometry')
%  @endverbatim


% defaults
prec = 8; % real numbers precision
path = '.';
misc.assigndefaults(varargin{:});

% Rename the information from fou struct
nreal = fou.nreal;
ndim = fou.ndim;
xyz_r = fou.xyz_r;
xyz_i = fou.xyz_i;
vecprod_r = fou.vecprod_r;
vecprod_i = fou.vecprod_i;
nm1 = fou.nm1;
nm2 = fou.nm2;

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
% modes for xyz
fwrite(fid,4, 'int');
fwrite(fid,nm1, 'int');
fwrite(fid,4, 'int');
% modes for vecprod
fwrite(fid,4, 'int');
fwrite(fid,nm2, 'int');
fwrite(fid,4, 'int');
% real part lagrangian mesh
for imode = 1:nm1
for i=1:ndim; 
   fwrite(fid,8*nreal, 'int');
   fwrite(fid,xyz_r(:,i,imode),   'double');
   fwrite(fid,8*nreal, 'int');
end
end
% imag part lagrangian mesh
for imode = 1:nm1
for i=1:ndim; 
   fwrite(fid,8*nreal, 'int');
   fwrite(fid,xyz_i(:,i,imode),   'double');
   fwrite(fid,8*nreal, 'int');
end
end
% real part vector product 
for imode = 1:nm2
for i=1:ndim; 
   fwrite(fid,8*nreal, 'int');
   fwrite(fid,vecprod_r(:,i,imode),   'double');
   fwrite(fid,8*nreal, 'int');
end
end
% imag part vector product
for imode = 1:nm2
for i=1:ndim; 
   fwrite(fid,8*nreal, 'int');
   fwrite(fid,vecprod_i(:,i,imode),   'double');
   fwrite(fid,8*nreal, 'int');
end
end

fclose(fid);

return 
end

