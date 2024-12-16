function [mesh,sub_mesh] = Mesh_GetSubset(mesh,vol_data,seg_tag,d2_thresh,it)
% Create a new index in the mesh input structures for the mesh faces that
% are within a d2_thresh from the vol_data of the segmentation of given tag

if nargin < 5
    it = 1;
end

% Work with f,v
v = mesh.v;
f = mesh.f;

% Get edge of the volume subset by subtracting original volume with and
% eroded (by 1 pixel) volume
vol_subset = (vol_data.img(:,:,:,it)==seg_tag);
[xk,yk,zk] = ndgrid(-1:1);
se_cross = strel(sqrt(xk.^2 + yk.^2 + zk.^2) <=1);

vol_eroded = imerode(vol_subset,se_cross);
vol_subset = (vol_subset & ~vol_eroded);

% Get coordinates of volume subset
X_subset = vol_data.X(vol_subset);
Y_subset = vol_data.Y(vol_subset);
Z_subset = vol_data.Z(vol_subset);

% Save these points for later
sub_mesh.px = X_subset;
sub_mesh.py = Y_subset;
sub_mesh.pz = Z_subset;

% Debug to check subset makes sense
% checkX = [min(X_subset(:)) max(X_subset(:)); min(v(:,1)) max(v(:,1))]
% checkY = [min(Y_subset(:)) max(Y_subset(:)); min(v(:,2)) max(v(:,2))]
% checkZ = [min(Z_subset(:)) max(Z_subset(:)); min(v(:,3)) max(v(:,3))]

% Split in parts to avoid out of memory!
% Do not use more than 2 additional GB: 2*1024*1024*1024/8 ~ 0.25e9
max_nv = round(0.25e9/length(X_subset));
v_parts = ceil(size(v,1)/max_nv);
disp(['Splitting point clouds distance computation in ' ...
    num2str(v_parts) ' parts.']);
i_end   = [(1:(v_parts-1))*max_nv size(v,1)];
i_start = [1 i_end(1:end-1)+1];

dist2v = ones(size(v,1),1);
for i = 1:v_parts
    % For each vertex find the distance to each volume point
    d2temp = bsxfun(@minus,v(i_start(i):i_end(i),1),X_subset(:)').^2 ...
           + bsxfun(@minus,v(i_start(i):i_end(i),2),Y_subset(:)').^2 ...
           + bsxfun(@minus,v(i_start(i):i_end(i),3),Z_subset(:)').^2;
    
    % And keep only the minimum distance for each vertex
    dist2v(i_start(i):i_end(i)) = min(d2temp,[],2);
end

% Heavy smooth distance to keep clean edges (4 passes)
[conn,~,~] = meshconn(f,length(v));
dist2v = smoothsurf(dist2v,[],conn,4,1);

% Get value for faces
min_dist2 = median(dist2v(f),2);

% Check using d2_thresh (using d^2)
mesh.isinsubset = (min_dist2 <= d2_thresh);

% Create sub_set mesh
sub_mesh.f = f(mesh.isinsubset,:);
sub_mesh.d2 = min_dist2(mesh.isinsubset);
% [sub_mesh.v,sub_mesh.f] = removeisolatednode(v,sub_mesh.f);
sub_f = sub_mesh.f;

for i = 1:size(mesh.v_mov,3)
    % Do one mesh at a time
    [sub_mesh.v,sub_mesh.f] = removeisolatednode(mesh.v_mov(:,:,i),sub_f);
    % Close the lid with bad triangles...
    sub_mesh = Mesh_CloseSurface(sub_mesh);

%     WriteToVTU(['dbg_submesh_mov' num2str(i,'%02i') '.vtu'],...
%         sub_mesh.v,sub_mesh.f);
    sub_mesh.v_mov(:,:,i) = sub_mesh.v;
end

sub_mesh.nt = mesh.nt;

% WriteToVTU('dbg_LAA_v.vtu',v,f,'laa_i',mesh.isinsubset,'dist',min_dist2);
% 
% WriteToVTU('dbg_LAA_xyz.vtu',[X_subset, Y_subset, Z_subset],...
%     1:1:length(X_subset),'cell_type',ones([1 length(X_subset)]),...
%         1:1:length(X_subset))

