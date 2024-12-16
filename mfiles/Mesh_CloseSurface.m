function mesh = Mesh_CloseSurface(mesh)

f = mesh.f;
v = mesh.v;

if isfield(mesh,'face_type')
    % Make sure to use the same amount of faces actually available
    ft = mesh.face_type(1:size(f,1));
else
    ft = zeros(size(f,1),1);
end

% Obtain open edges (vertix indices not ordered)
openedge = surfedge(f);

% Find and separate edges
loops = extractloops(openedge);

i_end = find(isnan(loops))-1;
i_start = [1 i_end(1:end-1)+2];

n_loops = length(i_end);

for i_loop = 1:n_loops
    % Notice that loops include one point twice (at the beginning and at
    % the end of the array)
    this_loop = loops(i_start(i_loop):i_end(i_loop))';
    mesh.edges_v_ind{i_loop} = this_loop(1:end-1);
    lid_center = mean(v(this_loop(1:end-1),:),1);
   
    v = [v; lid_center];
    new_v = size(v,1);
    
    new_fs = [this_loop(2:end) this_loop(1:end-1)...
        repmat(new_v,size(this_loop(1:end-1)))];
        
    % Add the new triangular faces
    f = [f; new_fs];
    
    % Give these faces new indices
    ft = [ft; i_loop*ones(size(new_fs,1),1)];
    
end

mesh.v = v;
mesh.f = f;
mesh.face_type = ft;
