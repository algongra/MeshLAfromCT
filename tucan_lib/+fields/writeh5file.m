function writeh5file(fnm,varargin)

h5fid = H5F.create(fnm, 'H5F_ACC_TRUNC','H5P_DEFAULT','H5P_DEFAULT');
typ = 'H5T_NATIVE_DOUBLE';
for i = 1:2:length(varargin)

   varnm = varargin{i};
   var = varargin{i+1};

   h5dims = fliplr(size(var));
   rank = length(h5dims);
   sid = H5S.create_simple(rank , h5dims, h5dims);
   did = H5D.create(h5fid, varnm, typ, sid,'H5P_DEFAULT');
   H5D.write(did,'H5ML_DEFAULT','H5S_ALL','H5S_ALL','H5P_DEFAULT',var);
   H5S.close(sid);
   H5D.close(did);
end

H5F.close(h5fid)

return
end
