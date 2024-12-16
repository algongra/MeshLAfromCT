function overwriteh5dset(varlist,fnm)

h5fid = H5F.open(fnm, 'H5F_ACC_RDWR','H5P_DEFAULT');

typbyclass.double = 'H5T_NATIVE_DOUBLE';
typbyclass.single = 'H5T_NATIVE_FLOAT';
typbyclass.int32 = 'H5T_NATIVE_INT';
% ux

for i = 1:length(varlist)
   varname = varlist{i};
   var     = evalin('caller',varname);
   h5dims = fliplr(size(var));
   did = H5D.open(h5fid, varname,'H5P_DEFAULT');
   H5D.write(did,'H5ML_DEFAULT','H5S_ALL','H5S_ALL','H5P_DEFAULT',var);
   H5D.close(did);
end

H5F.close(h5fid)

return
end
