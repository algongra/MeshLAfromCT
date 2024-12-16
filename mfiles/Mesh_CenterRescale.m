function mesh = Mesh_CenterRescale(mesh,scaling_factor)
disp('Rescaling and centering mesh...');

% Check if center was computed before
if isfield(mesh,'center')
    % Reuse that variable
    center = mesh.center;
    disp('Using previously computed center');
else
    % center mesh points at point x=0 y=0 z=0:
    center = (max(max(mesh.v_mov,[],1),[],3)+min(min(mesh.v_mov,[],1),[],3))./2;
end

% Shift to have center at the origin
mesh.v_mov(:,1,:) = mesh.v_mov(:,1,:)-center(1);
mesh.v_mov(:,2,:) = mesh.v_mov(:,2,:)-center(2);
mesh.v_mov(:,3,:) = mesh.v_mov(:,3,:)-center(3);

% rescale distances and areas
mesh.v_mov = mesh.v_mov*scaling_factor;
mesh.areas_mov = mesh.areas_mov*scaling_factor.^2;

% Reassing first frame to v list
mesh.v = mesh.v_mov(:,:,1);

% Store translation/scaling parameters
mesh.center = center;
mesh.scaling_factor = scaling_factor;
end
