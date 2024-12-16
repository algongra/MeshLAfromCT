function mesh = Mesh_MoveWithCPD(temp_mesh,start_frame,nt_seg,folder_out)

% CPD Options.
opts.corresp=0;
opts.normalize=1;
opts.max_it=150;
opts.tol=1e-5;
opts.viz=0;
opts.outliers=0.1;
opts.fgt=2;
% Beta: Width of the Gaussian.
% Lambda: Higher = more regularization, but poorer fit.
opts.beta=2;
opts.lambda = 3;

% use nonrigid registration with lowrank kernel approximation
opts.method = 'nonrigid_lowrank';

% leave only 70 large eigenvectors/values to approximate G
opts.numeig = 70;

% use FGT to find the largest eigenvectore/values
opts.eigfgt = 2;

% Copy starting frame
mesh = temp_mesh(start_frame);

% Allocate moving vertices
dump_v_mov = zeros([size(mesh.v) nt_seg]);

mesh.nt = nt_seg;

% Store starting frame in new array v_mov
%%%% Prepare vars for parfor loop

% Moving start mesh
movin_mesh = mesh.v;

% Output variable
dump_v_mov(:,:,start_frame) = mesh.v;

parfor i = 1:nt_seg
    disp([' ################## Registering Frame ' num2str(i) ' ##################']);
    if i ~= start_frame
        % Perform the CPD registration:
        % displace the start_frame to match the other independent mesh
        [Transform, ~] = cpd_register(temp_mesh(i).v,...
            movin_mesh,opts);

        ParSave(fullfile(folder_out,['CPD_dump' num2str(i) '.mat']),Transform);
        
        dump_v_mov(:,:,i) = Transform.Y;
    end
end
mesh.v_mov = dump_v_mov;
folder_out_ori = folder_out;
save(fullfile(folder_out,'temp_registration.mat'),'mesh','folder_out_ori');
