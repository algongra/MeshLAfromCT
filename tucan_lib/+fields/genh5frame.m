function genh5frame(ux,uy,uz,p,phi,fnm)


h5fid = H5F.create(fnm, 'H5F_ACC_TRUNC','H5P_DEFAULT','H5P_DEFAULT');
typ = 'H5T_NATIVE_DOUBLE';
% ux
h5dims = fliplr(size(ux));
rank = length(h5dims);
sid = H5S.create_simple(rank , h5dims, h5dims);
did = H5D.create(h5fid, 'ux', typ, sid,'H5P_DEFAULT');
H5D.write(did,'H5ML_DEFAULT','H5S_ALL','H5S_ALL','H5P_DEFAULT',ux);
H5S.close(sid);
H5D.close(did);

% uy
h5dims = fliplr(size(uy));
rank = length(h5dims);
sid = H5S.create_simple(rank , h5dims, h5dims);
did = H5D.create(h5fid, 'uy', typ, sid,'H5P_DEFAULT');
H5D.write(did,'H5ML_DEFAULT','H5S_ALL','H5S_ALL','H5P_DEFAULT',uy);
H5S.close(sid);
H5D.close(did);

% uz
h5dims = fliplr(size(uz));
rank = length(h5dims);
sid = H5S.create_simple(rank , h5dims, h5dims);
did = H5D.create(h5fid, 'uz', typ, sid,'H5P_DEFAULT');
H5D.write(did,'H5ML_DEFAULT','H5S_ALL','H5S_ALL','H5P_DEFAULT',uz);
H5S.close(sid);
H5D.close(did);

% p
h5dims = fliplr(size(p));
rank = length(h5dims);
sid = H5S.create_simple(rank , h5dims, h5dims);
did = H5D.create(h5fid, 'p', typ, sid,'H5P_DEFAULT');
H5D.write(did,'H5ML_DEFAULT','H5S_ALL','H5S_ALL','H5P_DEFAULT',p);
H5S.close(sid);
H5D.close(did);

% phi
h5dims = fliplr(size(phi));
rank = length(h5dims);
sid = H5S.create_simple(rank , h5dims, h5dims);
did = H5D.create(h5fid, 'phi', typ, sid,'H5P_DEFAULT');
H5D.write(did,'H5ML_DEFAULT','H5S_ALL','H5S_ALL','H5P_DEFAULT',phi);
H5S.close(sid);
H5D.close(did);

H5F.close(h5fid)

return
end
