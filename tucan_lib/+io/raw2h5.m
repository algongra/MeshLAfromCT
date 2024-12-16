function raw2h5(basenm, ifr, varargin)
% @author M.Moriche
% @brief Function to transform one old frame to a
%        hdf5 format frame
%
% @date 29-12-2013 by M.Moriche \n
%       Documented
%
% @details
%
% NO OUTPUT
% 
% MANDATORY ARGUMENTS
% - basenm:
% - ifr:
%
% OPTIONAL ARGUMENTS
% - path: path where the old frame is saved (DEFAULT: './')
% - path2: path where the frame is to be saved (DEFAULT: './')
% - basenm2: basenm of the frame to be saved. (DEFAULT: basenm)
% - ifr2: frame index of the frame to be saved. (DEFAULTL: ifr)
% - items: Items in frame to be read and saved.
% - fnm2: file name of the outputframe. If set, overrides path2, basnm2...
%
% EXAMPLES:
%
% @code
% raw2h5('nacamoving',100) 
% raw2h5('nacamoving',100,'path','/data/IBcode') 
% raw2h5('nacamoving',100,'basenm2', 'ncmvng','ifr2',0,'path2','/otherdir') 
% @endcode

%

path='.';
items = {'ux','uy','p','phi','uxeast','uyeast'};
basenm2 = '*';
ifr2 = '*';
path2 = '*';
fnm2 = '*';

misc.assigndefaults(varargin{:});

if path2 == '*'
   path2 = path;
end
if ifr2 == '*'
   ifr2 = ifr;
end
if basenm2 == '*'
   basenm2 = basenm;
end

fr = io.leefullframe(basenm, ifr,'items',items,'path',path);

fr.ux = fr.ux(2:end-1,2:end-1);
fr.uy = fr.uy(2:end-1,2:end-1);
fr.p  = fr.p(2:end-1,2:end-1);
fr.phi  = fr.phi(2:end-1,2:end-1);
fr.uxeast = fr.uxeast(2:end-1);
fr.uyeast = fr.uyeast(2:end-1);

fnm = fullfile(path2, basenm2, sprintf('%s.%04d.h5frame',basenm2, ifr));

if ~strcmp(fnm2,'*')
   fnm = fnm2;
end

h5fid = H5F.create(fnm, 'H5F_ACC_TRUNC','H5P_DEFAULT','H5P_DEFAULT');


for i = 1:length(items) 
   fld = items{i};
   typ = 'H5T_NATIVE_DOUBLE';
   var = getfield(fr,fld);
   h5dims = fliplr(size(var));
   rank = length(h5dims);
   sid = H5S.create_simple(rank , h5dims, h5dims);
   did = H5D.create(h5fid, fld, typ, sid,'H5P_DEFAULT');
   H5D.write(did,'H5ML_DEFAULT','H5S_ALL','H5S_ALL','H5P_DEFAULT',var);
   H5S.close(sid);
   H5D.close(did);
end

H5F.close(h5fid)

return
end
