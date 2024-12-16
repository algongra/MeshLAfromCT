function parin = file2parin(fnm)
% @author A.Gonzalo
%
% @brief Function to read the binary file that contains the information
%        to impose a parabolic velocity profile in the inlet of the pulmonary
%        veins (PVs).
%
% @date 12-12-2019 by A.Gonzalo \n
%                  Created and documented
%
%
% MANDATORY ARGUMENTS
%  - fnm : filename (absolute path). [char].

fid=fopen(fnm,'r','l');
% read number of cork planes
dummy = fread(fid,1, 'int');
ncpls  = fread(fid,1, 'int');
dummy = fread(fid,1, 'int');
% read number of points in each vein
for i=1:4;
    dummy = fread(fid, 1, 'int');
    nvpts(i) = fread(fid,  1,'int');
    dummy = fread(fid, 1, 'int');
end
% read values of the velocity in each vein point normalized with the maximum
% velocity of each vein
for i=1:4; 
   dummy = fread(fid, 1, 'int');
   u_in{i} = fread(fid, nvpts(i), 'double');
   dummy = fread(fid, 1, 'int');
end

fclose(fid);
