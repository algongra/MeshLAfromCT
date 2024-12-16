function dim = getdimensions(fnm, varnm);

h5fid = H5F.open(fnm,'H5F_ACC_RDONLY','H5P_DEFAULT');

did = H5D.open(h5fid,varnm);
sid = H5D.get_space(did);
[rnk dim] = H5S.get_simple_extent_dims(sid);
dim = fliplr(dim);
H5S.close(sid);
H5D.close(did);
H5F.close(h5fid);

return
end
