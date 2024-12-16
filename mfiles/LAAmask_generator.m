function LAAarea = LAAmask_generator(folder_out,seg_tag,data,mesh,nmodes,nmodes1,i_meshing,T,nx,ny,nz,varargin)

% path_filename = [folder_out '/full_preproc.mat'];    
subfolder = 'LAAmesh';

%%%%% paths for additional routines
% CFD Solver matlab API
addpath('../tucan_lib');
% Meshing tools
addpath('../iso2mesh');
% Load NIfTI images 
addpath('../NIfTI_load');
% Folder with the CPD registration functions
addpath(genpath('../ucsd-cvil_tools'))


% defaults
seg_time = [];
assigndefaults(varargin{:});
ifcubicspln = false;
if ~isempty(seg_time)
    ifcubicspln = true;
    nt_seg = size(data.img,4); dt = T/nt_seg;
    seg_time = 0:dt:1-dt;
end

seg_index = i_meshing+1;

% Go back to the mm units and center of the original CT scan
% This modifies .v_mov (not .v!)
mesh = Mesh_UndoCenterRescale(mesh);

% Allow up to 2*dz distance
d2_thresh = (3*data.dz).^2;

% Work on the frame of interest
if ifcubicspln
   ttarget = seg_time(seg_index);
   tmesh = 0:mesh.dt:T-mesh.dt;
   % Finding times in tmesh closer to ttarget
   tdiff = abs(tmesh-ttarget);
   % If ttarget time value exist in tmesh: no need for interpolation
   if any(tdiff==0);
      it_mesh = find(tdiff==0);
      mesh.v = mesh.v_mov(:,:,it_mesh);
   % If ttarget time value does exist in tmesh: linear interpolation...
   % ACHTUNG!!! Do we need to use a higher order interpolation to obtain
   % a reasonably good LAA triangulated mesh?
   else
      tdiffsign = tmesh-ttarget;
      % time closer to ttarget (from right)
      il = find((tdiffsign + tdiff)==0, 1,'last'); tl = tmesh(il);
      % time closer to ttarget (from right)
      ir = find((tdiffsign - tdiff)==0, 1,'first'); tr = tmesh(ir);
      if ir-il ~= 1
        error('Corrupted seg_time. Array must have ascendent order');
      end
      % Create interpolation weights
      wl = (tr-ttarget)/(tr-tl);
      wr = (ttarget-tl)/(tr-tl);
      % Interpolate linearly
      mesh.v = wl.*mesh.v_mov(:,:,il) + wr.*mesh.v_mov(:,:,ir);
   end
else
   mesh.v = mesh.v_mov(:,:,seg_index);
end
if size(data.img,4) == 1
    % if read from the CT image, data.img has only the one frame needed
    it = 1; 
else
    % if loaded from matlab volume data, pick the corresponding time frame
    it = seg_index;
end
% Extract sub set mesh and add tag to original mesh (working on .v)
[mesh,LAAmesh] = Mesh_GetSubset(mesh,data,seg_tag,d2_thresh,it);

% WriteToVTU('temp0_LAA_closed.vtu',LAAmesh.v,LAAmesh.f);
% WriteToVTU('temp0_LAA_v.vtu',mesh.v,mesh.f);

% Bring back to original postion the old mesh
mesh = Mesh_CenterRescale(mesh,0.1);
mesh.v = mesh.v_mov(:,:,seg_index);

% as well as the new one
LAAmesh = Mesh_ShiftRescale(LAAmesh,mesh.center,0.1);
LAAmesh.px = 0.1*(LAAmesh.px-LAAmesh.center(1));
LAAmesh.py = 0.1*(LAAmesh.py-LAAmesh.center(2));
LAAmesh.pz = 0.1*(LAAmesh.pz-LAAmesh.center(3));

% Get Fourier coefficients (also areas, normals)
LAAmesh.dt = T/LAAmesh.nt;
LAAmesh = Mesh_SmoothFourier2(LAAmesh,nmodes);

mkdir([folder_out '/' subfolder]);

%%% DEBUG dump vtu files
% for it = 1:size(LAAmesh.v_mov,3)
%     LAAmesh.qual_mov(:,it) = meshquality(LAAmesh.v_mov(:,:,it),LAAmesh.f);
%     WriteToVTU([folder_out '/' subfolder '/LAA_TriMoving_'...
%         num2str(length(LAAmesh.f(:,1))) 'tri_' num2str(it,'%02i') '.vtu'],...
%         LAAmesh.v_mov(:,:,it),LAAmesh.f,...
%         'areas',LAAmesh.areas_mov(:,it),...
%         'quality',LAAmesh.qual_mov(:,it),...
%         'face_i',LAAmesh.face_type,...
%         'tri_i',1:size(LAAmesh.f,1),...
%         'normal_v',LAAmesh.normals_mov(:,1,it),...
%         LAAmesh.normals_mov(:,2,it),...
%         LAAmesh.normals_mov(:,3,it));
% end
% Save whole LAAmesh structure
save([folder_out '/' subfolder '/LAAmesh.mat'],'LAAmesh');

% Points used to find LAA
WriteToVTU([folder_out '/' subfolder '/LAA_xyz.vtu'],...
    [LAAmesh.px, LAAmesh.py, LAAmesh.pz],1:1:length(LAAmesh.px),...
    'cell_type',ones([1 length(LAAmesh.px)]),1:1:length(LAAmesh.px));

LAAarea = sum(LAAmesh.areas_mov(LAAmesh.face_type~=1,:));

% some things that we need to create a mask for postprocessing
caras = LAAmesh.f;
coeffs = LAAmesh.v_hat;
norm_out = LAAmesh.norm_out;
[~,sort_ind] = sort(cellfun(@length,LAAmesh.edges_v_ind));
edges_v_ind = LAAmesh.edges_v_ind(sort_ind);
edge_ind = edges_v_ind{end};
filename = fullfile([folder_out '/' ...
    sprintf('data_for_LAAmask_%ix%ix%i',nx,ny,nz)]);
save(filename,'caras','coeffs','nmodes1','norm_out','edge_ind')


