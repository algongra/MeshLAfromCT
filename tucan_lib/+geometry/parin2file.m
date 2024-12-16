function parin2file(ncpls,nvpts,u_in,fnm,varargin)
% @author A.Gonzalo
%
% @brief Function to generate the binary file that contains the information
%        to impose a parabolic velocity profile in the inlet of the pulmonary
%        veins (PVs).
%
% @date 12-12-2019 by A.Gonzalo \n
%                  Created and documented
%
%
% MANDATORY ARGUMENTS
%  - ncpls: number of planes used to create the corks of the veins. [integer].
%  - nvpts: number of points used to mesh each vein. [integer]. Size: 1xnvein,
%           where nvein is the numner of PVs.
%  - u_in: values of the velocity in each vein point normalized with the maximum
%          velocity of each vein to impose the parabolic profile in the veins.
%          [cell]. Size: 1xnvein, where nvein is the numner of PVs.
%  - fnm : filename. [char].

% OPTIONAL ARGUMENTS
%  - prec: precision of real numbers. [integer]. Default: prec = 8.
%  - path: path to save the file. [char]. Default: path = '.'.

% defaults
prec = 8;
path = '.';

misc.assigndefaults(varargin{:});

% sanity checks
if ~isinteger(ncpls)
   if isfloat(ncpls)
      ncpls = int32(ncpls);
   else
      error('ncpls must be an integer or a float');
      return
   end
elseif ~isinteger(nvpts)
   if isfloat(nvpts)
      nvpts = int32(nvpts);
   else
      error('nvpts must be an integer or a float');
      return
   end
elseif ~strcmp(class(u_in),'cell')
   error('u_in must be a cell');
   return
elseif ~strcmp(class(fnm),'char')
   error('fnm must be a char');
   return
elseif ~isinteger(prec) & ~isfloat(prec)
   error('prec must be an integer or a float');
   return
elseif ~strcmp(class(path),'char')
   error('path must be a char');
   return
elseif length(nvpts)~=length(u_in)
   error('nvpts and u_in must have the same length');
   return
end

% create path folder in case it does not exist yet
system(sprintf('mkdir -p %s',path));

% get full path of the filename
fullfnm = fullfile(path, fnm);

% open filename
fid=fopen(fullfnm,'w','l');
% write number of cork planes
fwrite(fid,4,    'int');
fwrite(fid,ncpls,'int');
fwrite(fid,4,    'int');
% write number of points in each vein
for iv=1:length(nvpts);
    fwrite(fid,4,        'int');
    fwrite(fid,nvpts(iv),'int');
    fwrite(fid,4,        'int');
end
% write values of the velocity in each vein point normalized with the maximum
% velocity of each vein
for iv=1:length(nvpts);
    fwrite(fid,prec*nvpts(iv),'int');
    fwrite(fid,u_in{iv},      'double');
    fwrite(fid,prec*nvpts(iv),'int');
end
% close filename
fclose(fid);

return 
end

