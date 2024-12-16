function open_bc = CT_FindOpenPlanes(data,tags,hstencil,i_meshing)

DEBUGdump_cloud = false;
if nargin < 3
    hstencil = 3;
end
if nargin < 4
    i_meshing = 1;
end

disp(sprintf('\nFinding open planes...'));
% Dilate LA
try
    se_cross = strel('sphere',3);
    se_ball2 = strel('sphere',5);
catch
    [xk,yk,zk] = ndgrid(-3:3);
    se_cross = strel(xk.^2 + yk.^2 + zk.^2<=3.^2);
    [xk,yk,zk] = ndgrid(-5:5);
    se_ball2 = strel(xk.^2 + yk.^2 + zk.^2<=5.^2);
end

LA_dil = imdilate(data.img(:,:,:,i_meshing)==tags.LA,se_cross);

% Loop for all PVs
n_pvs = length(tags.PVs);

%%%%%% START DEBUG
if DEBUGdump_cloud
    dumpPVcloud = zeros(size(LA_dil)); 
    dumpLAcloud = zeros(size(LA_dil));
end
%%%%%% END DEBUG

fprintf('Inlet.....');
for i_pv = 1:n_pvs
    fprintf('\b\b%02i',i_pv);
    % Some smoothing and filling needed
    if hstencil
        h3 = ones(hstencil*[1 1 1])./(hstencil.^3);
        PVmask = imfilter(data.img(:,:,:,i_meshing)==tags.PVs(i_pv),h3);
        PVmask = imfill(PVmask>0.5,'holes');
    else
        PVmask = data.img(:,:,:,i_meshing)==tags.PVs(i_pv);
    end

    % Dilate each PV
    PV_dil = imdilate(PVmask,se_cross);
    % PV_dil = imdilate(data.img(:,:,:,i_meshing)==tags.PVs(i_pv),se_cross);
    
    % Find intersecting volumes to find outward direction
    plane_voxelsPV = LA_dil & (data.img(:,:,:,i_meshing)==tags.PVs(i_pv));
    plane_voxelsLA = PV_dil & (data.img(:,:,:,i_meshing)==tags.LA);

    %    try
    %        voxPV = gpuArray(plane_voxelsPV);
    %        voxPV = imclose(voxPV,se_ball2);
    %        plane_voxelsPV = gather(vox_PV); 
    %        clear voxPV
    %
    %        voxLA = gpuArray(plane_voxelsLA);
    %        voxLA = imclose(voxLA,se_ball2); 
    %        plane_voxelsLA = gather(vox_LA); 
    %        clear voxLA
    %    catch
    %        plane_voxelsPV = imclose(plane_voxelsPV,se_ball2);
    %        plane_voxelsLA = imclose(plane_voxelsLA,se_ball2); 
    %    end

    %%%%%% START DEBUG
    if DEBUGdump_cloud
        dumpPVcloud = dumpPVcloud + i_pv*plane_voxelsPV;
        dumpLAcloud = dumpLAcloud + i_pv*plane_voxelsLA;
    end
    %%%%%% END DEBUG
    % Point in PV on the edge
    edge_cloud(i_pv).x = data.X(imdilate(plane_voxelsPV,se_cross));
    edge_cloud(i_pv).y = data.Y(imdilate(plane_voxelsPV,se_cross));
    edge_cloud(i_pv).z = data.Z(imdilate(plane_voxelsPV,se_cross));
    x_PV_edge = mean(data.X(plane_voxelsPV));
    y_PV_edge = mean(data.Y(plane_voxelsPV));
    z_PV_edge = mean(data.Z(plane_voxelsPV));
    % Point in LA on the edge
    x_LA_edge = mean(data.X(plane_voxelsLA));
    y_LA_edge = mean(data.Y(plane_voxelsLA));
    z_LA_edge = mean(data.Z(plane_voxelsLA));
    % Vector pointing outward (not exactly normal)
    v_out(1) = x_PV_edge-x_LA_edge;
    v_out(2) = y_PV_edge-y_LA_edge;
    v_out(3) = z_PV_edge-z_LA_edge;
    
    % Find intersecting volume
    plane_voxels = LA_dil & PV_dil;
    %{
    pl_PV{i_pv} = bwconncomp(plane_voxels);
    if pl_PV{i_pv}.NumObjects > 1
       [~,iStruct] = max(cellfun(@length,pl_PV{i_pv}.PixelIdxList));
       plane_voxels = 0*plane_voxels;
       plane_voxels(pl_PV{i_pv}.PixelIdxList{iStruct}) = 1;
       clear iStruct;
    end
    %}
    % Get voxel coordinates
    x_plane = data.X(plane_voxels);
    y_plane = data.Y(plane_voxels);
    z_plane = data.Z(plane_voxels);
    % Find plane center
    x_c(i_pv) = mean(x_plane);
    y_c(i_pv) = mean(y_plane);
    z_c(i_pv) = mean(z_plane);
    % Find Normal to plane by best fit
    norm_plane = affine_fit([x_plane y_plane z_plane]);
    % Choose orientation using v_out direction
    if v_out*norm_plane < 0
        norm_plane = - norm_plane;
    end
    nx_plane(i_pv) = norm_plane(1);
    ny_plane(i_pv) = norm_plane(2);
    nz_plane(i_pv) = norm_plane(3);
    
    max_r(i_pv) = sqrt( max( (x_c(i_pv)-x_plane).^2 +...
        (y_c(i_pv)-y_plane).^2 + (z_c(i_pv)-z_plane).^2));
    
end
% Repeat the same for the LV to get Mitral Valve
fprintf('\nMV...');
% Some LV smoothing and filling needed

if hstencil
    h3 = ones(hstencil*[1 1 1])./(hstencil.^3);
    LVmask = imfilter(data.img(:,:,:,i_meshing)==tags.LV,h3);
    LVmask = imfill(LVmask>0.5,'holes');
else
    LVmask = data.img(:,:,:,i_meshing)==tags.LV;
end

LV_dil = imdilate(LVmask,se_cross);

% Find intersecting volumes to find outward direction
plane_voxelsLV = LA_dil & (data.img(:,:,:,i_meshing)==tags.LV);
plane_voxelsLA = LV_dil & (data.img(:,:,:,i_meshing)==tags.LA);
% Image closing to avoid holes
%    try
%        voxLV = gpuArray(plane_voxelsLV);
%        voxLV = imclose(voxLV,se_ball2);
%        plane_voxelsLV = gather(vox_LV); 
%        clear voxLV
%
%        voxLA = gpuArray(plane_voxelsLA);
%        voxLA = imclose(voxLA,se_ball2); 
%        plane_voxelsLA = gather(vox_LA); 
%        clear voxLA
%    catch
        plane_voxelsPV = imclose(plane_voxelsPV,se_ball2);
        plane_voxelsLA = imclose(plane_voxelsLA,se_ball2); 
%    end
%%%%%% START DEBUG
if DEBUGdump_cloud
    dumpPVcloud = dumpPVcloud + (n_pvs+1)*plane_voxelsLV;
    dumpLAcloud = dumpLAcloud + (n_pvs+1)*plane_voxelsLA;
end
%%%%%% END DEBUG
% Point in PV on the edge
edge_cloud(n_pvs+1).x = data.X(imdilate(plane_voxelsLV,se_cross));
edge_cloud(n_pvs+1).y = data.Y(imdilate(plane_voxelsLV,se_cross));
edge_cloud(n_pvs+1).z = data.Z(imdilate(plane_voxelsLV,se_cross));
x_LV_edge = mean(data.X(plane_voxelsLV));
y_LV_edge = mean(data.Y(plane_voxelsLV));
z_LV_edge = mean(data.Z(plane_voxelsLV));
% Point in LA on the edge
x_LA_edge = mean(data.X(plane_voxelsLA));
y_LA_edge = mean(data.Y(plane_voxelsLA));
z_LA_edge = mean(data.Z(plane_voxelsLA));
% Vector pointing outward (not exactly normal)
v_out(1) = x_LV_edge-x_LA_edge;
v_out(2) = y_LV_edge-y_LA_edge;
v_out(3) = z_LV_edge-z_LA_edge;

% Find intersecting volume
plane_voxels = LA_dil & LV_dil;
% Get voxel coordinates
x_plane = data.X(plane_voxels);
y_plane = data.Y(plane_voxels);
z_plane = data.Z(plane_voxels);
% Find plane center
x_c(n_pvs+1) = mean(x_plane);
y_c(n_pvs+1) = mean(y_plane);
z_c(n_pvs+1) = mean(z_plane);
% Find Normal to plane by best fit
norm_plane = affine_fit([x_plane y_plane z_plane]);

% Choose orientation using v_out direction
if v_out*norm_plane < 0
    norm_plane = - norm_plane;
end

nx_plane(n_pvs+1) = norm_plane(1);
ny_plane(n_pvs+1) = norm_plane(2);
nz_plane(n_pvs+1) = norm_plane(3);

max_r(n_pvs+1) = sqrt( max( (x_c(n_pvs+1)-x_plane).^2 +...
        (y_c(n_pvs+1)-y_plane).^2 + (z_c(n_pvs+1)-z_plane).^2));

fprintf('\nDone.\n');

%%%%%% START DEBUG
if DEBUGdump_cloud
WriteToVTR(['../PointsCloud_' datestr(datetime('now'),'ddHHMMSS')...
    '.vtr'],data.x,data.y,data.z,'LAcloud',dumpLAcloud,'PVMVcloud',dumpPVcloud);
end
%%%%%% END DEBUG

% plot3(x_plane,y_plane,z_plane,'.'); hold on;
% quiver3(x_c(n_pvs+1),y_c(n_pvs+1),z_c(n_pvs+1),...
%     nx_plane(n_pvs+1),ny_plane(n_pvs+1),nz_plane(n_pvs+1),5,'k','LineWidth',3);
% 
% quiver3(x_c(n_pvs+1),y_c(n_pvs+1),z_c(n_pvs+1),...
%         v_out(1),v_out(2),v_out(3),5,'r','LineWidth',3);
%     
open_bc.x_c      = x_c     ;
open_bc.y_c      = y_c     ;
open_bc.z_c      = z_c     ;
open_bc.nx_plane = nx_plane;
open_bc.ny_plane = ny_plane;
open_bc.nz_plane = nz_plane;
open_bc.max_r    = max_r;
open_bc.edge_cloud    = edge_cloud;
