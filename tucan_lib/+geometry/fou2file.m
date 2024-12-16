function fou2file(xyz,vecprod,nm1,nm2,fnm,varargin)
% @author MGV
%
%
% MANDATORY ARGUMENTS
%  - xyz(n,ndim,nm1): matrix with Fourier modes  Lagrangian grid points [double]
%  - vecprod(n,ndim,nm2): matrix with Fourier modes of vector product  [double]
%  - nm1 : number of modes to be written for xyz
%  - nm2 : number of modes to be written for vecprod
%  - fnm : filename

% OPTIONAL ARGUMENTS
%  - prec: precision of real numbers. [integer] = 8
%  - path: path to save the file [string] = '.'
%
%

% defaults
prec = 8;          % real numbers precision
path = '.';

misc.assigndefaults(varargin{:});

[nreal ndim dummy] = size(xyz);

% I use the 2 to reconstruct using sines and cosines
xyz_r = 2*real(squeeze(xyz(:,:,1:nm1)));
xyz_i = -2*imag(squeeze(xyz(:,:,1:nm1)));

vecprod_r = 2*real(squeeze(vecprod(:,:,1:nm2)));
vecprod_i = -2*imag(squeeze(vecprod(:,:,1:nm2)));

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

