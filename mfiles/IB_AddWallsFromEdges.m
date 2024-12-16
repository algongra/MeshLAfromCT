function ib_mesh = IB_AddWallsFromEdges(ib_mesh,mesh,corklen,d_plane)
% Uses the mesh.v_edge_ordered field to build an additional set of points
% representing an inflow pipe for each PV

% Number of PVs (last edge is the MV)
n_pv = mesh.n_edges-1;

% Number of points that will be added
% n_add = sum(mesh.edges_typ<=n_pv);

% Allocate new arrays
% new_ib_pts = zeros(n_add,3,mesh.nt);
% new_ib_normals = zeros(n_add,3,mesh.nt);
% new_ib_areas = zeros(n_add,mesh.nt);
% new_ib_type = mesh.edges_typ(mesh.edges_typ<=n_pv);


% Work on all time frames
for it = 1:mesh.nt
    % Initialize empty
    new_pts = [];
    new_normals = [];
    new_areas = [];
    new_type = [];
    for i_pv = 1:n_pv
        % Identify triangle faces for current vein plane
        face_type_indices = (mesh.face_type == i_pv);
        
        % Take one normal using the forced normals from the ib_mesh
        normal_to_inlet = median(ib_mesh.normals(face_type_indices,:,it),1);
        
        % delta displacement
        dist_vec = repmat(d_plane,[1 3]);
        
        for iplane = 1:corklen
            % Displace all the first vertices of the edge
            new_pts    = [new_pts;...
                bsxfun(@plus,mesh.v_mov(mesh.v_edge_ordered{i_pv}(:,1),:,it),...
                iplane.*normal_to_inlet.*dist_vec)];
            
            % Estimate areas by taking length of each segment and
            % multiplying by length of dist_vec
            new_areas = [new_areas;...
                sqrt( sum( (mesh.v_mov(mesh.v_edge_ordered{i_pv}(:,1),:,it)-...
                mesh.v_mov(mesh.v_edge_ordered{i_pv}(:,2),:,it)).^2 ,2)).* ...
                d_plane];
            
            % Assign normal to inlet to these points too
            new_normals = [new_normals;...
                repmat(normal_to_inlet,...
                [size(mesh.v_edge_ordered{i_pv}(:,1)) 1])];
            
            % Assign type
            new_type = [new_type; ...
                (10+i_pv).*ones(size(mesh.v_edge_ordered{i_pv}(:,1)))];
                        
        end        
    end
    add_pts(:,:,it) = new_pts;
    add_normals(:,:,it) = new_normals;
    add_areas(:,it) = new_areas;
    % vecprod value obtained from area and normal
    add_vecprod(:,:,it) = 2.*bsxfun(@times,new_areas,new_normals);
end
% This does not change with time (store only once at the end)
add_type = new_type;

% Put everything back in ib_mesh
ib_mesh.pts = [ib_mesh.pts; add_pts];
ib_mesh.areas = [ib_mesh.areas; add_areas];
ib_mesh.normals = [ib_mesh.normals; add_normals];
ib_mesh.vecprod = [ib_mesh.vecprod; add_vecprod];
ib_mesh.type = [ib_mesh.type; add_type];

% If more than 1 frame, _hat are coeffs in Fourier space
if isfield(ib_mesh, 'pts_hat')
    % Compute Fourier transform and add velocity
    add_pts_hat = fft(add_pts,[],3)./size(add_pts,3);
    
    % Time
    nt  = mesh.nt;
    dt  = mesh.dt;
    
    % Wavenumber array
    k   = 2*pi*(-1/(2*dt):(1/(dt*nt)):(1/(2*dt)));
    kt  = ifftshift(k(1:end-1));
    kt  = repmat(reshape(kt,[1 1 mesh.nt]),[size(add_pts_hat,1) size(add_pts_hat,2) 1]);
    
    % Obtain velocities
    add_vel_hat = 1i.*kt.*add_pts_hat;
    add_vel = ifft(add_vel_hat,[],3);
    
    % additional FFTs
    add_areas_hat = fft(add_areas,[],2)./size(add_areas,2);
    add_normals_hat = fft(add_normals,[],3)./size(add_normals,3);
    add_vecprod_hat = fft(add_vecprod,[],3)./size(add_normals,3);
    
    % Put everything back in ib_mesh
    ib_mesh.pts_hat = [ib_mesh.pts_hat; add_pts_hat];
    ib_mesh.vel_hat = [ib_mesh.vel_hat; add_vel_hat];
    ib_mesh.areas_hat = [ib_mesh.areas_hat; add_areas_hat];
    ib_mesh.normals_hat = [ib_mesh.normals_hat; add_normals_hat];
    ib_mesh.vecprod_hat = [ib_mesh.vecprod_hat; add_vecprod_hat];
    ib_mesh.vel = [ib_mesh.vel; add_vel];
end

end
