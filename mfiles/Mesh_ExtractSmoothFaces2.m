function mesh = Mesh_ExtractSmoothFaces2(mesh,open_bc,faces_cyl_thresh,folder_out)

% Work better with short name local variable
v = mesh.v;
f = mesh.f;
normals = mesh.normals;

% faces_cyl_thresh = sqrt(mesh.mean_area)/4;

disp(['Locating in/outlets with faces_cyl_thres=' num2str(faces_cyl_thresh)]);
%%%%%%%%%%%%%%%%%%% Find Walls and Inlets/Outlets %%%%%%%%%%%%%%%%%%%%%%%%%
% Work on face centers instead of vertices!
% centers_x = (v(f(:,1),1)+v(f(:,2),1)+v(f(:,3),1))/3;
% centers_y = (v(f(:,1),2)+v(f(:,2),2)+v(f(:,3),2))/3;
% centers_z = (v(f(:,1),3)+v(f(:,2),3)+v(f(:,3),3))/3;

% number of open bc to process
n_bc = length(open_bc.x_c);

% Logical array for wall faces and each inlet/outlet
wall_face_all = true(size(f(:,1)));
inout_face_all(:,:) = false(size(f(:,1),1),n_bc);
face_type = zeros(size(f(:,1)));

% Compute distances in portions to avoid OUT OF MEMORY
max_nv = 1e4;
nv = size(v,1);
n_loops = floor(nv/max_nv);

% Loop among all open bc
for i = 1:n_bc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For each vertex of all triangles find the minimum distance from the
    % point cloud defining the inlet/outlet
    dist2v_cloud = zeros(nv,1);
    % Loop in partitions to avoid OUT OF MEMORY
    for il = 1:n_loops
        indices = ((il-1)*max_nv+1):(il*max_nv);
        dist2v = bsxfun(@minus,v(indices,1),open_bc.edge_cloud(i).x').^2 ...
               + bsxfun(@minus,v(indices,2),open_bc.edge_cloud(i).y').^2 ...
               + bsxfun(@minus,v(indices,3),open_bc.edge_cloud(i).z').^2;
        dist2v_cloud(indices) = min(dist2v,[],2);
    end
    
    % Last (or only until end of vector) until end 
    indices = ((n_loops)*max_nv+1):nv;
    dist2v = bsxfun(@minus,v(indices,1),open_bc.edge_cloud(i).x').^2;
    dist2v = dist2v + bsxfun(@minus,v(indices,2),open_bc.edge_cloud(i).y').^2;
    dist2v = dist2v + bsxfun(@minus,v(indices,3),open_bc.edge_cloud(i).z').^2;
    dist2v_cloud(indices) = min(dist2v,[],2);

    % Smooth distance from cloud of points using smoothsurf
    [conn,~,~] = meshconn(f,length(v));
    dist2v_cloud = smoothsurf(dist2v_cloud,[],conn,4,1);

    % Assign to each triangle the minimum distance among all their vertices
    % min_dist2 = min(dist2v_cloud(f),[],2);
    % Assign to each triangle the mean distance among all their vertices
    min_dist2 = median(dist2v_cloud(f),2);
    
    min_d2{i} = min_dist2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Vector from each triangle vertex to cylinder base center
    % First triangle vertex
    vec_x(:,1) = v(f(:,1),1)-open_bc.x_c(i);
    vec_y(:,1) = v(f(:,1),2)-open_bc.y_c(i);
    vec_z(:,1) = v(f(:,1),3)-open_bc.z_c(i);
    % Second triangle vertex
    vec_x(:,2) = v(f(:,2),1)-open_bc.x_c(i);
    vec_y(:,2) = v(f(:,2),2)-open_bc.y_c(i);
    vec_z(:,2) = v(f(:,2),3)-open_bc.z_c(i);
    % Third triangle vertex
    vec_x(:,3) = v(f(:,3),1)-open_bc.x_c(i);
    vec_y(:,3) = v(f(:,3),2)-open_bc.y_c(i);
    vec_z(:,3) = v(f(:,3),3)-open_bc.z_c(i);
        
    % Normalize normal vector defining cylinder orientation
    nor_mag = sqrt(open_bc.nx_plane(i).^2+open_bc.ny_plane(i).^2+open_bc.nz_plane(i).^2);
    nor_x = open_bc.nx_plane(i)./nor_mag;
    nor_y = open_bc.ny_plane(i)./nor_mag;
    nor_z = open_bc.nz_plane(i)./nor_mag;
    
    % Find the distance along the normal for each vertex, then pick the
    % maximum
    dot_prod = vec_x.*nor_x + vec_y.*nor_y + vec_z.*nor_z;
    dp{i} = abs(min(dot_prod,[],2));
    % dot_prod = (vec_x.*repmat(nor_x,[1 3]))+...
    %     (vec_y.*repmat(nor_y,[1 3]))+(vec_z.*repmat(nor_z,[1 3]));
    
    
    % Find the distance on the plane
    radius = sqrt((vec_x-dot_prod.*nor_x).^2+...
        (vec_y-dot_prod.*nor_y).^2+...
        (vec_z-dot_prod.*nor_z).^2);
    
    face_dot_prod = open_bc.nx_plane(i)*normals(:,1) + ...
        open_bc.ny_plane(i)*normals(:,2) + ...
        open_bc.nz_plane(i)*normals(:,3);
    face_dp{i} = face_dot_prod;
    % Choose faces wihtin a thresh thick cylinder with radius the max distance
    % from the centers and that have similar normal orientation
    % (|dot|<0.5 => angle<pi/6)
    %if i == n_bc
    if 1
        % Treatment for MV
        keep_face = min_dist2 >= faces_cyl_thresh; % |...            
            % (abs(min(dot_prod,[],2))>=sqrt(faces_cyl_thresh));
            % (face_dot_prod<0.3) |...
    else
        % Treatment for PVs
        keep_face = min_dist2 >= faces_cyl_thresh;
    end

    %%%%%%%%%%%%% Remove disconnected faces belonging to edge %%%%%%%%%%%%%
    % Temporary array for in_out faces
    f_io_temp       = reshape(f(repmat(~keep_face,[1 3])),[],3);
    discon_faces    = finddisconnsurf(f_io_temp);    
    flip_faces = [];        
    if length(discon_faces) > 1
        disp([num2str(length(discon_faces))...
            ' disconnected faces found (' num2str(i) ') !']);
        % Reorder subset (descending)
        [~,I] = sort(cellfun(@length,discon_faces),'descend');
        discon_faces = discon_faces(I);
        % Assume the first(largest) subset is the correct one
        for i_d = 2:length(discon_faces)
            % Loop among all the first vertices present in the wrong face list
            for i_dd = 1:size(discon_faces{i_d},1)
                % One vertex (1st) is enough to identify out-of-place faces
                flip_faces = [flip_faces;...
                    find(any(f==discon_faces{i_d}(i_dd,1),2))];
            end
        end
        % Change the disconnected faces that contain these vertices
        % from not_keep(wall) to keep(inout)
        keep_face(flip_faces) = true;
    end    
    
    %%%%%%%%%%%%% Remove disconnected faces belonging to wall %%%%%%%%%%%%%
    % Temporary array for wall faces
    f_wa_temp       = reshape(f(repmat(keep_face,[1 3])),[],3);
    discon_faces    = finddisconnsurf(f_wa_temp);
    flip_faces = [];    
    % Remove disconnected faces
    if length(discon_faces) > 1
        disp([num2str(length(discon_faces))...
            ' disconnected faces found (' num2str(i) ') !']);
        % Reorder subset (descending)
        [~,I] = sort(cellfun(@length,discon_faces),'descend');
        discon_faces = discon_faces(I);
        % Assume the first(largest) subset is the correct one
        for i_d = 2:length(discon_faces)
            % Loop among all the first vertices present in the wrong face list
            for i_dd = 1:size(discon_faces{i_d},1)
                % One vertex (1st) is enough to identify out-of-place faces
                flip_faces = [flip_faces;...
                    find(any(f==discon_faces{i_d}(i_dd,1),2))];
            end
        end
        % Change the disconnected faces that contain these vertices
        % from not_keep(wall) to keep(inout)
        keep_face(flip_faces) = false;
    end

    keep_face = Mesh_RemoveLooseFaces(f,keep_face);
    
    wall_face_all = wall_face_all & keep_face;
    
    inout_face_all(:,i) = ~keep_face;
    
    face_type = face_type + i*(~keep_face);
end
%end

%  % Put Wall and Inout faces in different arrays
%  f_wall = f(repmat(wall_face_all,[1 3]));
%  f_wall = reshape(f_wall,[],3);
%  f_inout_all = f(repmat(~wall_face_all,[1 3]));
%  f_inout_all = reshape(f_inout_all,[],3);
%  
%  % Single inlet/outlet for each cell
%  for i = 1:n_bc
%      f_inout{i} = f(repmat(inout_face_all(:,i),[1 3]));
%      f_inout{i} = reshape(f_inout{i},[],3);
%  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
WriteToVTU([folder_out 'LAmesh_full_prefix_' datestr(datetime('now'),'ddHHMMSS')...
    '.vtu'],v,f,'face_i',face_type,...
    'dist_1',min_d2{1},'face_dp1',face_dp{1},'dp1',dp{1},...
    'dist_2',min_d2{2},'face_dp2',face_dp{2},'dp2',dp{2},...
    'dist_3',min_d2{3},'face_dp3',face_dp{3},'dp3',dp{3},...
    'dist_4',min_d2{4},'face_dp4',face_dp{4},'dp4',dp{4},...
    'dist_5',min_d2{5},'face_dp5',face_dp{5},'dp5',dp{5});
%}
% Base inputs
flnm = [folder_out 'LAmesh_full_prefix_'...
        datestr(datetime('now'), 'ddHHMMSS') '.vtu'];
base_args = {flnm, v, f, 'face_i', face_type};

% Initialize a cell array to hold all arguments
args = base_args;

% Dynamically add inputs based on the length of min_d2
for i = 1:length(min_d2)
    args = [args, {['dist_' num2str(i)], min_d2{i}, ...
                   ['face_dp' num2str(i)], face_dp{i}, ...
                   ['dp' num2str(i)], dp{i}}];
end

% Call the function with dynamically constructed arguments
WriteToVTU(args{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nChecking open edges...\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Additional Check fo MV %%%%%%%%%%%%%%%%%%%%%%%
fix_acc = [];
fix_count = 0;
redo_edges = true;
%type_mvpv = 5;
type_mvpv = max(face_type);
while redo_edges
    f_mv = f(face_type==type_mvpv,:);

    % Obtain open edges (vertix indices not ordered)
    openedge_all = surfedge(f_mv);

    % Find and separate edges
    loops = extractloops(openedge_all);

    % Each edge has a NaN at the end: split them
    i_end = find(isnan(loops))-1;
    i_start = [1 i_end(1:end-1)+2];

    % Check the length of the edges: if less than 10 vertices, drop the edge
    i_diff  = i_end-i_start;
    i_sta_fix = i_start(i_diff<10);
    i_end_fix = i_end(i_diff<10);
    n_fix = length(i_sta_fix);

    n_edges = length(i_end);
    if n_edges>1 && fix_count<10
        fix_count = fix_count + 1;
        warning('Locating Inlets/Outles having problems. Trying to fix it!')
        for i = 1:n_fix
            fix_loop = unique(loops(i_sta_fix(i):i_end_fix(i)));
            fix_f_mask = sum(sum(bsxfun(@eq,sort(f,2),reshape(fix_loop,[1 1 length(fix_loop)])),3),2)==3;
            fix_f = find(fix_f_mask);
            ft_before = face_type(fix_f)
            
            % If == , changes it to max of the neighbours
            % else becomes zero
            if face_type(fix_f) == type_mvpv
                face_type(fix_f) = 0;
            else
                face_type(fix_f) = type_mvpv;
            end
            
            ft_after = face_type(fix_f)
            fix_acc = [fix_acc; squeeze(fix_f)];
        end
    else
        redo_edges = false;
    end
end
%%%%%%%%%%%%%%%%%%%%% End of additional Check of MV %%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
WriteToVTU([folder_out 'LAmesh_full_prefix_AA_' datestr(datetime('now'),'ddHHMMSS')...
    '.vtu'],v,f,'face_i',face_type,...
    'dist_1',min_d2{1},'face_dp1',face_dp{1},'dp1',dp{1},...
    'dist_2',min_d2{2},'face_dp2',face_dp{2},'dp2',dp{2},...
    'dist_3',min_d2{3},'face_dp3',face_dp{3},'dp3',dp{3},...
    'dist_4',min_d2{4},'face_dp4',face_dp{4},'dp4',dp{4},...
    'dist_5',min_d2{5},'face_dp5',face_dp{5},'dp5',dp{5});
%}
% Base inputs
flnm = [folder_out 'LAmesh_full_prefix_AA'...
        datestr(datetime('now'), 'ddHHMMSS') '.vtu'];
base_args = {flnm, v, f, 'face_i', face_type};

% Initialize a cell array to hold all arguments
args = base_args;

% Dynamically add inputs based on the length of min_d2
for i = 1:length(min_d2)
    args = [args, {['dist_' num2str(i)], min_d2{i}, ...
                   ['face_dp' num2str(i)], face_dp{i}, ...
                   ['dp' num2str(i)], dp{i}}];
end

% Call the function with dynamically constructed arguments
WriteToVTU(args{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





fix_count = 0;
redo_edges = true;
while redo_edges
    f_wall = f(face_type==0,:);

    % Obtain open edges (vertix indices not ordered)
    openedge_all = surfedge(f_wall);

    % Find and separate edges
    loops = extractloops(openedge_all);

    % Each edge has a NaN at the end: split them
    i_end = find(isnan(loops))-1;
    i_start = [1 i_end(1:end-1)+2];

    % Check the length of the edges: if less than 10 vertices, drop the edge
    i_diff  = i_end-i_start;
    i_sta_fix = i_start(i_diff<10);
    i_end_fix = i_end(i_diff<10);
    n_fix = length(i_sta_fix);

    n_edges = length(i_end);
    if n_edges~=n_bc && fix_count<10
        fix_count = fix_count + 1;
        warning('Locating Inlets/Outles having problems. Trying to fix it!')
        for i = 1:n_fix
            fix_loop = unique(loops(i_sta_fix(i):i_end_fix(i)));
            fix_f_mask = sum(sum(bsxfun(@eq,sort(f,2),reshape(fix_loop,[1 1 length(fix_loop)])),3),2)==3;
            fix_f = find(fix_f_mask);
            ft_before = face_type(fix_f)
            
            % If ==0, changes it to max of the neighbours
            % else becomes zero
            face_type(fix_f) = (face_type(fix_f)==0)*...
                max(face_type(max(1,fix_f-5):min(fix_f+5,length(face_type))));
            
            ft_after = face_type(fix_f)
            fix_acc = [fix_acc; squeeze(fix_f)];
            disp(['Face(s) ' num2str(squeeze(fix_f)') ' changed from ' ...
                num2str(squeeze(ft_before)') ' to ' num2str(squeeze(ft_after)')]);
        end
    else
        redo_edges = false;
    end
end
if ~isempty(fix_acc)
    disp([num2str(length(fix_acc)) ' single faces fixed!']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{
    WriteToVTU([folder_out 'LAmesh_full_noproc_' datestr(datetime('now'),'ddHHMMSS')...
        '.vtu'],v,f,'face_i',face_type,...
        'dist_1',min_d2{1},'face_dp1',face_dp{1},'dp1',dp{1},...
        'dist_2',min_d2{2},'face_dp2',face_dp{2},'dp2',dp{2},...
        'dist_3',min_d2{3},'face_dp3',face_dp{3},'dp3',dp{3},...
        'dist_4',min_d2{4},'face_dp4',face_dp{4},'dp4',dp{4},...
        'dist_5',min_d2{5},'face_dp5',face_dp{5},'dp5',dp{5});
    disp('#############  OJO: SE HA GENERADO UN LAmesh_full_noproc!!!')
    %}
    % Base inputs
    flnm = [folder_out 'LAmesh_full_noproc_'...
            datestr(datetime('now'), 'ddHHMMSS') '.vtu'];
    base_args = {flnm, v, f, 'face_i', face_type};
    
    % Initialize a cell array to hold all arguments
    args = base_args;
    
    % Dynamically add inputs based on the length of min_d2
    for i = 1:length(min_d2)
        args = [args, {['dist_' num2str(i)], min_d2{i}, ...
                       ['face_dp' num2str(i)], face_dp{i}, ...
                       ['dp' num2str(i)], dp{i}}];
    end
    
    % Call the function with dynamically constructed arguments
    WriteToVTU(args{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


% figure;
fprintf('\nObtaining open edges......');
edges_pts = [];
edges_typ = [];
for i_edge = 1:n_edges
    fprintf('\b\b\b\b\b%02i/%02i',i_edge,n_edges);
    % i_edge
    % Vertices of current edge
    v_circle{i_edge} = loops(i_start(i_edge):i_end(i_edge));
    % Find the triangle edges connecting the points in the current edge
    logi = any(repmat(openedge_all(:,1),[1 size(v_circle{i_edge},2)])==...
        repmat(v_circle{i_edge},[size(openedge_all(:,1),1) 1]),2);
    edge{i_edge} = openedge_all(logi,:);
    % AGG_BEG
    % Avoid getting two edge segments when one vertex is shared by segments
    % from different edges
    % (may occur when different edges are too close [same triangle])
    %
    % Find vertices' indices of edge
    edge_v = unique(edge{i_edge}(:));
    % Count how many times each vertex's index appeares in the edge
    cnt_v_rep = sum(repmat(edge{i_edge}(:),[1 length(edge_v)])==...
                    repmat(edge_v',[length(edge{i_edge}(:)) 1]),1);
    % Create new edge discarding segments containing vertex's index with less
    % than two appearances
    if any(cnt_v_rep<2)
       edge{i_edge} = edge{i_edge}(~sum(ismember(edge{i_edge},...
                                                 edge_v(cnt_v_rep<2)),2),:);
    end
    % AGG_END
    
    % Order the points
    v_edge_ordered{i_edge} = orderloopedge(edge{i_edge});
    
    % To find the edge index, check the index of the triangles containing
    % one of the point in the edge (the first one...)
    % i_edge does not correspond to the cylinder index used before !!!
    edge_index_temp = face_type(any(f==v_edge_ordered{i_edge}(1,1),2));
    edge_index = mean(edge_index_temp(edge_index_temp~=0));
    
    % Put points coordinates into temporary arrays
    xv = v(v_edge_ordered{i_edge}(:,1),1);
    yv = v(v_edge_ordered{i_edge}(:,1),2);
    zv = v(v_edge_ordered{i_edge}(:,1),3);
    
    % Low pass filter to smooth close to circular edges
    % OLD, this created very wiggling curves
    % n_harm_keep = max([5 ceil(length(edge{i_edge})/10)]);
    n_harm_keep = 10;
    
    xv_h = fft(xv);
    yv_h = fft(yv);
    zv_h = fft(zv);
    
    xv_h(n_harm_keep+1:end-(n_harm_keep-1)) = 0;
    yv_h(n_harm_keep+1:end-(n_harm_keep-1)) = 0;
    zv_h(n_harm_keep+1:end-(n_harm_keep-1)) = 0;
    
    xv_s{i_edge} = ifft(xv_h);
    yv_s{i_edge} = ifft(yv_h);
    zv_s{i_edge} = ifft(zv_h);
    
    % Put smoothed points back into the original array of vertices
    v(v_edge_ordered{i_edge}(:,1),1) = xv_s{i_edge};
    v(v_edge_ordered{i_edge}(:,1),2) = yv_s{i_edge};
    v(v_edge_ordered{i_edge}(:,1),3) = zv_s{i_edge};
    
    edges_pts = [edges_pts; [xv_s{i_edge} yv_s{i_edge} zv_s{i_edge}] ];
    edges_typ = [edges_typ; edge_index*ones(size(xv_s{i_edge}))];
end
fprintf('\nDone.\n'); 

% Reorder v_edge_ordered using the face_type order
[~,edges_ordered] = sort(unique(edges_typ,'stable'));

% Put everything in a structure as function output argument
mesh.f = f;
mesh.face_type = face_type;
mesh.v = v;
mesh.edges_pts = edges_pts;
mesh.edges_typ = edges_typ;
mesh.n_edges = n_edges;
mesh.v_edge_ordered = v_edge_ordered(edges_ordered);
