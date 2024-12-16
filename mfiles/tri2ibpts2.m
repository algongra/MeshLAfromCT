function [mesh,ib_mesh] = tri2ibpts2(mesh,varargin)
% Goes from tringular mesh to centers of triangles, for each of the edge
% surfaces one normal direction is obtained with a best fit plane

%% Inlcude paths of external libraries
addpath('../tucan_lib');

%% Assign varargin or default values
% Set flag to specify if it is a moving mesh (ifmov=1) or a static one (ifmov=0)
ifmov = 1;
misc.assigndefaults(varargin{:});

dbg_on = false;

% Work on local variables copied from input structure
n_edges = mesh.n_edges;
f = mesh.f;
face_type = mesh.face_type;

if isfield(mesh,'v_mov')
    v_mov = mesh.v_mov;
    vel   = mesh.vel;
    nt = size(v_mov,3);
end

if ifmov == 0
    v_mov = mesh.v;
    vel   = zeros(size(v_mov));
    nt = mesh.nt;
end

% Remove old field if present
if isfield(mesh,'ib_pts')
    mesh = rmfield(mesh,'ib_pts');
end

if isfield(mesh,'ib_normals')
    mesh = rmfield(mesh,'ib_normals');
end

if isfield(mesh,'ib_vecprod')
    mesh = rmfield(mesh,'ib_vecprod');
end

if isfield(mesh,'ib_vel')
    mesh = rmfield(mesh,'ib_vel');
end
    

for it = 1:nt
    disp(sprintf('\nWorking on frame %i',it));
    % Work on normals for each frame separately
    if nt == 1
       % Only one frame (no moving normals)
       normals = mesh.normals;
	    vecprod = mesh.vecprod;
	    areas = mesh.areas;
    elseif isfield(mesh,'normals_mov')
        % More than one frame
        normals = mesh.normals_mov(:,:,it);
        vecprod = mesh.vecprod_mov(:,:,it);
        areas = mesh.areas_mov(:,it);
    else
        error('tri2ibpts2: normals_mov was not saved in the mesh structure.');
        return
    end
    
    
    % Normals to the each inlet/outlets set to be the same (best fit)
    normals_inout = zeros([n_edges 3]);
    for i_edge = 1:n_edges
        % Points centroid of each face (edge)
        edge_face_pts = meshcentroid(v_mov(:,:,it),f(face_type==i_edge,:));
        
        % Fit plane and obtain normal vector
        [normals_inout(i_edge,:),~,~] = affine_fit(edge_face_pts);
        
        % Compute other normal as median of all traingles' normals
        normals_med = median(normals(face_type==i_edge,:));
        
        % This is to check that we got the right normal from the plane fit: if
        % it's in the opposite direction (dot negative) just reverse it
        if dot(normals_med,normals_inout(i_edge,:))<0
            normals_inout(i_edge,:) = -normals_inout(i_edge,:);
        end
        % Save also to the full normals list
        normals(face_type==i_edge,1) = normals_inout(i_edge,1);
        normals(face_type==i_edge,2) = normals_inout(i_edge,2);
        normals(face_type==i_edge,3) = normals_inout(i_edge,3);

        % Save normals for cork present (20+i_edge)
        normals(face_type==20+i_edge,1) = normals_inout(i_edge,1);
        normals(face_type==20+i_edge,2) = normals_inout(i_edge,2);
        normals(face_type==20+i_edge,3) = normals_inout(i_edge,3);

        % Save vecprod from best-fit normal and actual areas
        vecprod(face_type==i_edge,1) = normals_inout(i_edge,1)*2*areas(face_type==i_edge);
        vecprod(face_type==i_edge,2) = normals_inout(i_edge,2)*2*areas(face_type==i_edge);
        vecprod(face_type==i_edge,3) = normals_inout(i_edge,3)*2*areas(face_type==i_edge);

        vecprod(face_type==20+i_edge,1) = normals_inout(i_edge,1)*2*areas(face_type==20+i_edge);
        vecprod(face_type==20+i_edge,2) = normals_inout(i_edge,2)*2*areas(face_type==20+i_edge);
        vecprod(face_type==20+i_edge,3) = normals_inout(i_edge,3)*2*areas(face_type==20+i_edge);
    end
    
    if dbg_on
%       save(['LAmesh_full_normals_' num2str(length(f(:,1))) ...
%           'tri_'  datestr(datetime('now'),'yyyymmdd_HHMM') '.mat'],...
%           'mesh');
       WriteToVTU(['LAmesh_full_normals_' num2str(length(f(:,1))) ...
           'tri_'  datestr(datetime('now'),'yyyymmdd_HHMM') '.vtu'],...
           mesh.v,mesh.f,'Areas',areas,...
           'mesh_qual',mesh.qual,'face_i',face_type,...
           'normal_v',normals(:,1),normals(:,2),normals(:,3));
    end
    
    % Obtain triangles' centers
    ib_pts = meshcentroid(v_mov(:,:,it),f);
    
    % Obtain triangle velocity
    ib_vel = (vel(f(:,1),:,it)+vel(f(:,2),:,it)+vel(f(:,3),:,it))/3;
    
    mesh.ib_pts(:,:,it) = ib_pts;
    mesh.ib_normals(:,:,it) = normals;
    mesh.ib_vecprod(:,:,it) = vecprod;
    mesh.ib_vel(:,:,it) = ib_vel;
end

% Copy to ib_* variables parameters which were not changed
if nt == 1
    % Only one frame (no changing areas/pts)
    mesh.ib_areas = mesh.areas;
    mesh.ib_type = face_type;
    % mesh.ib_vecprod = mesh.vecprod;    
elseif isfield(mesh,'areas_mov') && isfield(mesh,'areas_hat')
    
    % More than one frame: areas and point are changing
    mesh.ib_areas = mesh.areas_mov;
    mesh.ib_areas_hat = mesh.areas_hat;
    mesh.ib_type = face_type;

    % mesh.ib_vecprod = mesh.vecprod_mov;    

    % Compute Fourier Transform of pts and normals
    mesh.ib_pts_hat = fft(mesh.ib_pts,[],3)./size(mesh.ib_pts,3);
    mesh.ib_normals_hat = fft(mesh.ib_normals,[],3)./size(mesh.ib_normals,3);
    mesh.ib_vel_hat = fft(mesh.ib_vel,[],3)./size(mesh.ib_vel,3);
    mesh.ib_vecprod_hat = fft(mesh.ib_vecprod,[],3)./size(mesh.ib_vecprod,3);
else
    error('tri2ibpts2: areas_mov/hat was not saved in the mesh structure.');
    return
end


ib_mesh.pts     = mesh.ib_pts;
ib_mesh.vel     = mesh.ib_vel;
ib_mesh.areas   = mesh.ib_areas;
ib_mesh.normals = mesh.ib_normals;
ib_mesh.vecprod = mesh.ib_vecprod;

ib_mesh.type    = mesh.ib_type;

if isfield(ib_mesh, 'edges_pts')
    ib_mesh.edges_pts = mesh.edges_pts;
end
if isfield(ib_mesh, 'edges_typ')
    ib_mesh.edges_typ = mesh.edges_typ;
end


% If multiple frames are processes save also Fourier space data
if nt>1
    ib_mesh.pts_hat     = mesh.ib_pts_hat;
    ib_mesh.vel_hat     = mesh.ib_vel_hat;
    ib_mesh.areas_hat   = mesh.ib_areas_hat;
    ib_mesh.normals_hat = mesh.ib_normals_hat;
    ib_mesh.vecprod_hat = mesh.ib_vecprod_hat;
end
