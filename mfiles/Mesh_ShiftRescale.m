function mesh = Mesh_ShiftRescale(mesh,shift_xyz,scaling_factor)
disp('Rescaling and shifting mesh...');

if ~isfield(mesh,'v_mov')
    mesh.v_mov = mesh.v;
end

% Shift to have center at the origin
mesh.v_mov(:,1,:) = mesh.v_mov(:,1,:)-shift_xyz(1);
mesh.v_mov(:,2,:) = mesh.v_mov(:,2,:)-shift_xyz(2);
mesh.v_mov(:,3,:) = mesh.v_mov(:,3,:)-shift_xyz(3);

% rescale distances and areas
mesh.v_mov = mesh.v_mov*scaling_factor;
if isfield(mesh,'areas_mov')
    mesh.areas_mov = mesh.areas_mov*scaling_factor.^2;
end

% Reassing first frame to v list
mesh.v = mesh.v_mov(:,:,1);

% Store translation parameters
mesh.center = shift_xyz;
mesh.scaling_factor = scaling_factor;
end
