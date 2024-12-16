function mesh = mask2mesh(data,i_meshing,mesh_keepratio,smooth_iter,MaxArea,tags,faces_cyl_thresh, folder_out)
% Extract the mesh, from a level set function. The CT segmentation is still
% used to identify the inlets and outlet

%% User parameter
% Decide how many elements to keep
% mesh_keepratio = 0.1;

% Segmentation path
%seg_path = '~/Dropbox/fluid-dynamics/ucsd_siemens/262829/seg-nii-lorenzo/';

% Frame used to obtain reference mesh (max LA volume around 40% RR)
% i_meshing = 8;

if nargin < 6
    % Inclusion tags
    tags.LA = 4;
    tags.LAA = 10;
    % Neighbouring tags
    tags.LV = 1;
    tags.PVs = 11:14;
end
if nargin < 7
    % Multyplying factor of mean_area of mesh triangles
    faces_cyl_thresh = 0.5;
end

% Find boundary planes
open_bc = CT_FindOpenPlanes(data,tags,0,i_meshing);

% Get the LA mask
data.mask = (data.img(:,:,:,i_meshing) == tags.LA ...
    | data.img(:,:,:,i_meshing) == tags.LAA);

% Obtain mesh (vertices and faces) and put in a working structure mesh
% 0: No additional smoothing on the volume (already applied in main script)
mesh = CT_MeshAndResample(data,1,mesh_keepratio,smooth_iter,0);

% Next steps are going to move the points defining the edges of each
% face: keep original also a copy (mesh_orig)
% mesh_orig = mesh;

% Make sure the orientation is correct, and add needed flag
% Now inside Mesh_AnalyzePM
% mesh = Mesh_CheckNormalsOrientation(mesh);

% Update faces areas, quality and normals
mesh = Mesh_AnalyzePM(mesh);

% Extract different faces from each one of the open bc defining
% in/outlets and smmoth the edges
% mesh = Mesh_ExtractSmoothFaces(mesh,open_bc,sqrt(mesh.mean_area));
disp(['Mean Area: ' num2str(mesh.mean_area)]);
disp([' setting treshold as ' num2str(faces_cyl_thresh) ' of mean area.']);
faces_cyl_thresh = faces_cyl_thresh*mesh.mean_area;
mesh = Mesh_ExtractSmoothFaces2(mesh,open_bc,faces_cyl_thresh,folder_out);

% Get faces areas, quality and normals
mesh = Mesh_AnalyzePM(mesh);
if nargin<6
    % Set MaxArea cutoff (used to limit IB points "volume")
    MaxArea = mean(mesh.areas)+2*std(mesh.areas);
end

% Split big triangles in 3 (not ideal: decreases mesh quality...)
mesh = Mesh_SetMaxArea(mesh,MaxArea);

% Update faces areas, quality and normals
mesh = Mesh_AnalyzePM(mesh);

%%% DEBUG
WriteToVTU([folder_out sprintf('debug_mesh.%02i.vtu',i_meshing)],...
           mesh.v,mesh.f,'Areas',mesh.areas,'mesh_qual',mesh.qual,...
           'face_i',mesh.face_type,'tri_i',1:size(mesh.f,1),...
           'normal_v',mesh.normals(:,1),mesh.normals(:,2),mesh.normals(:,3));
 
% 
% stlwrite([seg_path '../tri_mesh' num2str(i_meshing) '.stl'],mesh.f,mesh.v);
% 
% save([seg_path '../tri_mesh' num2str(i_meshing) '.mat'],'mesh');
return
% p = patch('Faces',mesh.f,'Vertices',mesh.v);
% p.FaceColor = 'red'; p.FaceAlpha = 0.5;
%%
