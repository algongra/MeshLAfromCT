function writeh5file(varlist,fnm)

h5fid = H5F.create(fnm, 'H5F_ACC_TRUNC','H5P_DEFAULT','H5P_DEFAULT');

typbyclass.single = 'H5T_NATIVE_FLOAT';
typbyclass.double = 'H5T_NATIVE_DOUBLE';
typbyclass.int32 = 'H5T_NATIVE_INT';
% ux

for i = 1:length(varlist)
   varname = varlist{i};
   var     = evalin('caller',varname);
   h5dims = fliplr(size(var));
   rank = length(h5dims);

   % modification for one dimensional arrays
   if length(var(:)) == size(var,1)
      rank = 1;
      h5dims = size(var,1);
   end

   classname = class(var);
   typ = typbyclass.(classname); 
   sid = H5S.create_simple(rank , h5dims, h5dims);
   did = H5D.create(h5fid, varname, typ, sid,'H5P_DEFAULT');
   H5D.write(did,'H5ML_DEFAULT','H5S_ALL','H5S_ALL','H5P_DEFAULT',var);
   H5S.close(sid);
   H5D.close(did);
end

H5F.close(h5fid)

return
end
