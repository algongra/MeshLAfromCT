function grid2file(xyz,vol,fnm,varargin)
% @author M.Moriche
%
% @brief Function to write coordinates of Lagrangian mesh defined by
%        the matrix xy and associated volume vol in file fnm.
%
% @date 3-12-2013 by M.Moriche \n
%       Created. Extracted from previous scripts.
%
% @details
%
% NO OUTPUT
%
% MANDATORY ARGUMENTS
%  - xy: matrix with coordinates of Lagrangian grid points [double]
%  - vol: associated volume of each Lagrangian marker. [double]
%  - fnm: file name to write. [string]
%
% OPTIONAL ARGUMENTS
%  - prec: precision of real numbers. [integer] = 8
%  - xcontrol: control point. Used for moving aspects in code 
%              subroutines. [double] = [0 0]
%  - path: path to save the file [string] = '.'
%
%  EXAMPLES:
%
%  @verbatim
%  grid2file(xy,vol,fnm)
%  grid2file(xy,vol,fnm,'xcontrol',[-1 0])
%  grid2file(xy,vol,fnm,'path','/data/geometry','xcontrol',[1 1])
%  @endverbatim
%

% defaults
prec = 8;          % real numbers precision
path = '.';

[nreal ndim] = size(xyz);
xcontrol = zeros(1,ndim);  % control point

misc.assigndefaults(varargin{:});

xcontrol,

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
% control point
for i=1:ndim 
fwrite(fid,8,     'int');
fwrite(fid,xcontrol(i),'double');
fwrite(fid,8,     'int');
end 
% lagrangian size
fwrite(fid,4, 'int');
fwrite(fid,nreal, 'int');
fwrite(fid,4, 'int');
nreal
% lagrangian mesh
for i=1:ndim; 
fwrite(fid,8*nreal, 'int');
fwrite(fid,xyz(:,i),   'double');
fwrite(fid,8*nreal, 'int');
end
% marker volume
fwrite(fid,8,   'int');
fwrite(fid,vol, 'double');
fwrite(fid,8,   'int');
%
fclose(fid);

return 
end

