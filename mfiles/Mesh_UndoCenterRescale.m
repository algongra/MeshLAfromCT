function mesh = Mesh_UndoCenterRescale(mesh)
disp('Restoring original scaling and position of mesh.v_mov...');

if isfield(mesh,'scaling_factor')
    scaling_factor = mesh.scaling_factor;
else
    scaling_factor = 0.1;
end

% rescale distances and areas
mesh.v_mov = mesh.v_mov/scaling_factor;
mesh.areas_mov = mesh.areas_mov/scaling_factor.^2;

% Shift to have center at the origin
mesh.v_mov(:,1,:) = mesh.v_mov(:,1,:)+mesh.center(1);
mesh.v_mov(:,2,:) = mesh.v_mov(:,2,:)+mesh.center(2);
mesh.v_mov(:,3,:) = mesh.v_mov(:,3,:)+mesh.center(3);

end
