function mesh = Mesh_AddCorks(mesh,corklen,d_plane,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Generation of corks (inside walls) %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input required:
% mesh: a structure that contains f,v_mov,...
% corklen: number of planes to be added for each cork
% d_plane: distance between each additional plane

%% Inlcude paths of external libraries
addpath('../tucan_lib');

%% Assign varargin or default values
% Set flag to specify if it is a moving mesh (ifmov=1) or a static one (ifmov=0)
ifmov = 1;
i_start_vein = 1;
misc.assigndefaults(varargin{:});

if nargin < 4
    i_start_vein = 1;
end

disp(sprintf('\nAdding corks...'));

if ifmov == 0
    nt = 1;
    mesh.v_mov = mesh.v;
    mesh.vel = zeros(size(mesh.v));
    mesh.qual_mov = mesh.qual;
    mesh.areas_mov = mesh.areas;
    mesh.normals_mov = mesh.normals;
    mesh.nt = nt;
end

if isfield(mesh,'v_mov')
    nt = size(mesh.v_mov,3);
end

% Getting number of veins
%   - Each vein is numbered with a different type
%   - The mitral valve is numbered with another type
nvein = mesh.n_edges-1;

for it = 1:nt
    % Work on 1 frame at a time
    mesh.v = mesh.v_mov(:,:,it);
    mesh.qual = mesh.qual_mov(:,it);
    mesh.areas = mesh.areas_mov(:,it);
    mesh.normals = mesh.normals_mov(:,:,it);
    
    % add points of corks located at veins inlets to triangular mesh:
    for i=i_start_vein:i_start_vein+nvein-1
        % Identify triangle faces for current vein plane
        face_type_indices = (mesh.face_type == i);
        inlet_faces = mesh.f(face_type_indices,:);
        
        % These quantities will be just copied
        new_areas = mesh.areas(face_type_indices,:);
        new_qual = mesh.qual(face_type_indices,:);
        new_normals = mesh.normals(face_type_indices,:);
        
        % Take unique sorted reference indices
        inlet_faces_uni = unique(inlet_faces);
        % new_vertices_ids = size(mesh.v,1)+(1:length(inlet_faces_uni));
        
        for iplane = 1:corklen
            % Map old vertices onto newly added vertices
            new_faces = inlet_faces(:);
            for ivert = 1:length(inlet_faces_uni)
                new_faces(new_faces==inlet_faces_uni(ivert)) = size(mesh.v,1)+ivert;
            end
            new_faces = reshape(new_faces,[length(new_faces)/3 , 3]);
            
            % Take one normal (median)
            normal_to_inlet = median(mesh.normals(face_type_indices,:),1);
            
            % Still using mm while dxs are in cm
            dist_vec = repmat(d_plane,[1 3]);
            % Displace all new vertices
            new_vertices    = bsxfun(@plus,mesh.v(inlet_faces_uni,:),iplane.*normal_to_inlet.*dist_vec);
            
            % Put everything back into structure
            % Copy these
            mesh.areas      = [mesh.areas; new_areas];
            mesh.qual      = [mesh.qual; new_qual];
            mesh.normals    = [mesh.normals; new_normals];
            mesh.v          = [mesh.v; new_vertices];
            
            if it==1
                % This entries are only added once during the time loop
                % These are newly computed (corks have indices of 100 + vein index)            
                mesh.f          = [mesh.f; new_faces];            
                mesh.face_type  = [mesh.face_type; (20+i)*ones(size(new_faces,1),1)];
            end
        end
    end
    % Work on 1 frame at a time
    mesh.v_new(:,:,it) = mesh.v;
    mesh.qual_new(:,it) = mesh.qual;
    mesh.areas_new(:,it) = mesh.areas;
    mesh.normals_new(:,:,it) = mesh.normals ;
end

mesh.v_mov = mesh.v_new;
mesh.qual_mov = mesh.qual_new;
mesh.areas_mov = mesh.areas_new;
mesh.normals_mov = mesh.normals_new;

mesh = rmfield(mesh,{'v_new','qual_new','areas_new','normals_new'});

end
