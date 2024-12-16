function body = file2body(fnm,varargin)
% @author M.Moriche
%
% @brief 
% @date 26-11-2014 by M.Moriche \n
%       Created from body2file
% @date 24-02-2014 by M.Moriche \n
%       Edited to obtain all the data
% @date 06-06-2016 by A.Gonzalo \n
%       Added normflag as optional argument
%
% @details
%
% OUTPUT
%  struct with fields of the lagrangian mesh
%
% MANDATORY ARGUMENTS
%  - fnm: file name to read. [string]
%
% OPTIONAL ARGUMENTS
%  - normflag: normal vectors flag. [logical]
%              (change it to true if normal vectors are included in fnm)
%
%  EXAMPLES:
%
%  @verbatim
%  bb = file2body(fnm)
%  bb = file2body(fnm,'normflag','true')
%  @endverbatim

% defaults
normflag = false;

misc.assigndefaults(varargin{:});

%
fid=fopen(fnm,'r','l');
% precision
dummy = fread(fid,1, 'int');
prec  = fread(fid,1, 'int');
dummy = fread(fid,1, 'int');
% dimensions
dummy = fread(fid, 1, 'int');
ndim = fread(fid,  1,'int');
dummy = fread(fid, 1, 'int');
% lagrangian size
dummy = fread(fid, 1, 'int');
nreal = fread(fid, 1, 'int');
dummy = fread(fid, 1, 'int');
% lagrangian mesh
xyz = zeros(nreal,ndim);

for i=1:ndim; 
   dummy = fread(fid, 1, 'int');
   xyz(:,i) = fread(fid, nreal, 'double');
   dummy = fread(fid, 1, 'int');
end
% marker volume
dummy = fread(fid,1,   'int');
vol = fread(fid,nreal, 'double');
dummy = fread(fid,1,   'int');
% normal vectors
if normflag
   normvec = zeros(nreal,ndim);

   for i=1:ndim;
       dummy = fread(fid, 1, 'int');
       normvec(:,i) = fread(fid, nreal, 'double');
       dummy = fread(fid, 1, 'int');
   end
end
% number of bodies
dummy = fread(fid, 1, 'int');
nbody = fread(fid, 1, 'int');
dummy = fread(fid, 1, 'int');
% control point
xyzc = zeros(nbody,ndim);

for i=1:ndim 
   fread(fid, 1,     'int');
   xyzc(:,i) = fread(fid, nbody,'double');
   fread(fid, 1,     'int');
end 
% begin index
fread(fid,1,'int');
ib = fread(fid,nbody,'int');
fread(fid,1,'int');
% end index
fread(fid,1,'int');
ie = fread(fid,nbody,'int');
fread(fid,1,'int');
% specific inertia (I/M)
ntensor = 0.5*ndim*(ndim+1);

fread(fid,1,'int');
spi = fread(fid,nbody*ntensor,'double');
fread(fid,1,'int');
% density ratio
dummy = fread(fid,1,    'int');
rhoratio = fread(fid,1, 'double');
dummy = fread(fid,1,    'int');
% body volume
dummy = fread(fid,1,    'int');
bodyvol = fread(fid,1, 'double');
dummy = fread(fid,1,    'int');


fclose(fid);

body.prec       = prec;
body.ndim       = ndim;
body.nreal      = nreal;
body.xyz        = xyz;
body.vol        = vol;
if normflag
   body.normvec = normvec;
end
body.xyzc       = xyzc;
body.ib         = ib;
body.ie         = ie;
body.rhoratio   = rhoratio;
body.bodyvol    = bodyvol;

return 
end

