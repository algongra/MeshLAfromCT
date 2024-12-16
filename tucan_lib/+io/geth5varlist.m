function varlist = geth5varlist(fnm)

varlist = {};

plist = 'H5P_DEFAULT';
fid = H5F.open(fnm,'H5F_ACC_RDONLY',plist);
gid = H5G.open(fid,'/');
root_info = H5G.get_info(gid);
idx_type = 'H5_INDEX_NAME';
order = 'H5_ITER_DEC';
for j = 0:root_info.nlinks-1
   obj_id = H5O.open_by_idx(fid,'/',idx_type,order,j,plist);
   obj_info = H5O.get_info(obj_id);
   nm = H5I.get_name(obj_id);
   switch(obj_info.type)
       case H5ML.get_constant_value('H5G_LINK')
           fprintf('Object #%d %s is a link.\n',j,nm);
       case H5ML.get_constant_value('H5G_GROUP')
           fprintf('Object #%d %s is a group.\n',j,nm);
       case H5ML.get_constant_value('H5G_DATASET')
           %fprintf('Object #%d %s is a dataset.\n',j,nm);
           varlist = [varlist strrep(nm,'/','')];
       case H5ML.get_constant_value('H5G_TYPE')
           fprintf('Object #%d %s is a named datatype.\n',j,nm);
   end
   H5O.close(obj_id);
end
H5G.close(gid);
H5F.close(fid);


return
end
