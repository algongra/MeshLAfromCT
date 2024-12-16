function mesh = Mesh_SmoothFourier2(mesh,n_harm_keep)

% Work on private variables
f = mesh.f;

disp('Performing FFTs...');

% FFT in time
v_hat = fft(mesh.v_mov,[],3);

% If 2 arguments drop some harmonics
if nargin > 1
    if n_harm_keep < size(v_hat,3)
        % Remove some higher frequencies
        v_hat(:,:,round(n_harm_keep)+1:end-(round(n_harm_keep)-1)) = 0;
        % Bring back to time domain
        mesh.v_mov = ifft(v_hat,[],3);
    end
end

% Save Fourier Coefficient in mesh structure as well
mesh.v_hat = v_hat./size(mesh.v_mov,3);

% Time
nt  = mesh.nt;
dt  = mesh.dt;
% FFT testing
% t   = 0:dt:dt*(nt-1);
% y   = sin(2*pi*t/(dt*nt));
% ype = 2*pi/(dt*nt)*cos(2*pi*t/(dt*nt));

% Wavenumber array
k   = 2*pi*(-1/(2*dt):(1/(dt*nt)):(1/(2*dt)));
kt  = ifftshift(k(1:end-1));
kt  = repmat(reshape(kt,[1 1 mesh.nt]),[size(v_hat,1) size(v_hat,2) 1]);

mesh.vel_hat = 1i.*kt.*v_hat;
mesh.vel = ifft(mesh.vel_hat,[],3);


% Allocate areas and normals in time
mesh.areas_mov = zeros(size(f,1),mesh.nt);
mesh.normals_mov = zeros(size(f,1),3,mesh.nt);
mesh.vecprod_mov = zeros(size(f,1),3,mesh.nt);

% Get mesh properties (area,quality and normals)
mesh = Mesh_AnalyzePM(mesh);

% Fourier Transform of areas and normals
mesh.areas_hat = fft(mesh.areas_mov,[],2)./size(mesh.areas_mov,2);
mesh.normals_hat = fft(mesh.normals_mov,[],3)./size(mesh.normals_mov,3);
mesh.vecprod_hat = fft(mesh.vecprod_mov,[],3)./size(mesh.vecprod_mov,3);
return
