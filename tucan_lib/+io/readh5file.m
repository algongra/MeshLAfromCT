function [fr, failitems] = readh5file(fnm, items, varargin)
% @author M.Moriche
% @date Created 06-06-13 by M.Moriche
% @date Modified 07-06-13 by M.Moriche
% @date 07-06-13 by M.Moriche
%       added extra arguments in H5D.read
% @date 22-07-13 by M.Moriche
%       added error handling and failitems as output
% @date 13-12-13 by M.Moriche
%       Organized and documented for trunk
% @date 25-12-13 by M.Moriche
%       sizebyvar as optional argument
% 
% @brief Function to read data from a file saved with HDF5 lib
%
% @details
%
% Uses low level function in HDF5 library. High level libraries are
%  not available in all Matlab releases.
%
% Mandatory arguments
%  - fnm: HDF5 file to be read.
%  - items: Datasets of the file to be read.
%
% Optional arguments.
%  - sizebyvar: struct where each field is a variable name
%               whose key is an struct with the following 
%               (MANDATORY) keys: see HDF5 MANUALS FOR REF.
%               -offset
%               -stride
%               -block
%               -count
% Examples:
%
%  @code
%  fr = readh5file('/data/nacamoving_102.0000.h5frame',{'ux','uy','xux'})
%  @endcode
%
%
%

sizebyvar = struct;
misc.assigndefaults(varargin{:});

failitems = {};

h5fid = H5F.open(fnm,'H5F_ACC_RDONLY','H5P_DEFAULT');

fr = struct;

for i = 1:length(items)
   varnm = items{i};
   try
      did = H5D.open(h5fid,varnm);
      sid = 'H5S_ALL';
      mid = 'H5S_ALL';
      uu = fields(sizebyvar);
      if find(strcmp(varnm , uu))
         ndim = length(sizebyvar.(varnm).offset);
         sid = H5D.get_space(did);
         H5S.select_hyperslab (sid, 'H5S_SELECT_SET', ...
                               fliplr(sizebyvar.(varnm).offset),...
                               fliplr(sizebyvar.(varnm).stride),...
                               fliplr(sizebyvar.(varnm).count),...
                               fliplr(sizebyvar.(varnm).block));
         mid = H5S.create_simple(ndim, fliplr(sizebyvar.(varnm).count), []);
      end

      fr.(items{i}) = H5D.read(did,...
       'H5ML_DEFAULT',mid, sid,'H5P_DEFAULT');
      H5D.close(did);
      if strcmp(class(fr.(items{i})),'char')
         fr.(items{i}) = deblank(fr.(items{i}));
      end

      if find(strcmp(varnm , uu))
         H5S.close(sid);
         H5S.close(mid);
      end
   catch
      failitems = [failitems items{i}];
   end
end
H5F.close(h5fid);
end
