clear; close; clc;

dbstop if error

%%%%% paths for additional routines
% CFD Solver matlab API
addpath('../tucan3v6_mfiles');
% Meshing tools
addpath('../iso2mesh');
% Load NIfTI images 
addpath('../NIfTI_load');
% Folder with the CPD registration functions
addpath(genpath('../ucsd-cvil_tools'));

timer_start=tic; disp('Starting now');

%% Full script from segmentation to TUCAN fortran files
% Name for case, will be used to save/load files
caseid = '27292';

% Extension of the segmentation files
seg_ext = '.nii.gz';

% times when segmentations were obtained
seg_time = [18 26 34 43 51 59 67 75 84 93]/100;
% dt use for spline interpolation
dt_spln = 0.025;

% Base folder where output folders will be created
% Klone
base_path = '/gscratch/mambrino/algongra/Projects/Hearts/ct_processing_la/';
%{
% TheGreat
base_path = '/Users/alex/Projects/Hearts/ct_processing_la_Klone/';
%}

% Folder containing the seg.nii.gz files with the segmentation
% Klone
seg_path = ['/gscratch/mambrino/algongra/Projects/Hearts/Bolus/'...
            'CT_scans_data/' caseid '/Cine_data/seg-nii/'];
%{
% TheGreat
seg_path = ['/Users/alex/Projects/Hearts/Bolus/'...
            'CT_scans_data/' caseid '/Cine_data/seg-nii/'];
%}
%{
seg_path = ['/Volumes/data28/alex/Projects/Hearts/Bolus/'...
            'CT_scans_data/' caseid '/Cine_data/seg-nii/'];
%}
disp(['Using segmentation in ' seg_path]);

% ---- Resolution

res = 256;

% ---- Patient specific tags (different than usual)

% Inclusion tags
tags.LA = 2;
tags.LAA = 3;
% Neighbouring tags
tags.LV = 1;
tags.PVs = 4:7;

% Cardiac period [s]
T = 1;
% Number of pulmunary veins
nvein = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% CFD solver parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting dx, dy and dz (dh) in cm!
%
nx =  res;  ny = nx;  nz = nx;
x0 = -6.5;  y0 = x0;  z0 = x0;
xf =  6.5;  yf = xf;  zf = xf;
Lc = xf-x0;
%
if nx == 144; resnm = 'lowres'; elseif nx == 256; resnm = 'nomres'; end;
%
dh = (xf - x0)/nx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Mesh USER PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ---- Segmentation volume smoothing parameters

% 4D smoothing on volumetric segmentation
smoothing_stencil = [5,5,5,1];
filter_type = 'gaussian';


% ---- Frame meshing parameters

% Decide how many elements of the initial mesh to keep
% +--------------------+-----------------------+----------------+
% | Spatial Resolution | Points per box length | mesh_keepratio |
% +--------------------+-----------------------+----------------+
% |        Low         |       nx = 144        |     0.1-0.3    |
% |      Nominal       |       nx = 256        |     0.3-0.5    |
% |        High        |       nx = 384        |     0.9-1.0    |
% +--------------------+-----------------------+----------------+
mesh_keepratio = 0.5;

% Number of smoothing iterations on triangular mesh (with iso2mesh)
smooth_iter = 2;

% Max Area of faces [mm^2] 
MaxArea = 1.5*(10*dh).^2;

% Multiplying factor of mean_area of mesh triangles as threshold
% to detect open faces of PV and MV
% [need to decrease with increasing resolution: 384 ~ 1, 256 ~ 0.5]
% Default value: 0.5
faces_cyl_thresh0 = 0.5;


% ---- Corks and fourier modes

% Number of planes to replicate to form the corks of the PVs [12]
corklen = 12;

% Distance between each cork plane [in mm like the whole CT/mesh processing]
d_plane = 10*dh;

% number of Fourier modes to be kept
nmodes = 10;
nmodes1 = nmodes-1;
nmodes2 = 2*nmodes1;

% ---- Output folders

% Folder name to save
casename = sprintf('%s_%imodes',caseid,nmodes);
folder_out = sprintf('%s/%s_%s_%g_%s/',base_path,casename,resnm,...
                     mesh_keepratio,datestr(datetime('now'),'yyyymmdd_HHMM'));

% Folder to save debug files
folder_debug = sprintf('%s/%s_%s_%g_debug/',base_path,casename,resnm,...
                       mesh_keepratio);

% ---- Restart files

% Info:
%  - reload_volume_path:       After segmentation volume smoothing is performed
%  - reload_mesh_path:         After meshing is performed
%  - reload_registration_path: After CPD registration is performed
%  - reload_corks_path:        After pulmonary veins corks are created
reload_volume_path = sprintf('%s/data_levset.mat',folder_debug);
reload_mesh_path = sprintf('%s/temp_mesh.mat',folder_debug);
reload_registration_path = sprintf('%s/temp_registration.mat',folder_debug);
reload_corks_path = sprintf('%s/temp_corks.mat',folder_debug);

% Flag to include double LA wall (two layers of points on the LA separated
% a distance dh in the normal direction)
ifdoublewall = true;



% Sanity check
if exist('faces_cyl_thresh_mod','var') & exist('i_mod','var') &...
   exist('ied_mod','var')
   if ~isequal(size(faces_cyl_thresh_mod),size(i_mod),size(ied_mod))
      % Reshape 2D matrices to 1D arrays
      faces_cyl_thresh_mod = reshape(faces_cyl_thresh_mod,...
                                     [1 numel(faces_cyl_thresh_mod)]);
      i_mod = reshape(i_mod,[1 numel(i_mod)]);
      ied_mod = reshape(ied_mod,[1 numel(ied_mod)]);
   else
      errtxt = ['faces_cyl_thresh_mod, i_mod, ied_mod arrays MUST have '...
                'the same size'];
      error(errtxt);
   end
end

% Create folder_out and folder_debug only if user is not using restart files
if ~exist('reload_volume_path') & ~exist('reload_mesh_path') &...
   ~exist('reload_registration_path') & ~exist('reload_corks_path')
   % Create output folder
   mkdir(folder_out);
   % Create debug folder
   mkdir(folder_debug);
end

%% -------------------  START PROCESSING -----------------

% Count number of frames to be read
nt_seg = length(dir([seg_path '*' seg_ext]));
disp([num2str(nt_seg) ' frames found!']);

% ---- Registrarion process parameters

%  Get initial segmentation frame to use
%  (obtained automatically as the file where LA volume is larger)
seg_tags = [tags.LA tags.LAA];
i_meshing = seg2vol_LA_fun(seg_path,seg_ext,seg_tags,nt_seg,T);

%% - - - - - - ITK VOLUME SMOOTHING - - - - - - -
tic
disp(' ----------- SMOOTHING ITK VOLUMES -----------')
% Smoothing segmentations
if exist('reload_volume_path')
   if isfile(reload_volume_path)
      disp(['Loading volume data from ' reload_volume_path]);
      load(reload_volume_path);
   end
else
   data = CT_GetSmoothSegWithLevSet(seg_path,nt_seg,seg_ext,tags,...
                                    smoothing_stencil,filter_type);

   % Save data for restart
   vars_to_save = {'caseid','seg_ext','seg_time','dt_spln','base_path',...
                   'seg_path','res','tags','T','nvein','res','nx','ny','nz',...
                   'x0','y0','z0','xf','yf','zf','resnm','dh',...
                   'smoothing_stencil','filter_type','mesh_keepratio',...
                   'smooth_iter','MaxArea','corklen','d_plane','nmodes',...
                   'nmodes1','nmodes2','casename','folder_out',...
                   'folder_debug','ifdoublewall'};
   if exist('faces_cyl_thresh_mod','var') & exist('i_mod','var') &...
      exist('ied_mod','var') 
      vars_to_save = [vars_to_save,'faces_cyl_thresh_mod','i_mod','ied_mod'];
   end
   vars_to_save = [vars_to_save,'nt_seg','seg_tags','i_meshing','data'];
   save(fullfile(folder_debug,'data_levset.mat'),vars_to_save{:},'-v7.3');
end

time_smoothing = toc;

%% - - - - - - MESHING INDIVIDUAL FRAMES - - - - - - - -

disp(' ----------- MESHING INDIVIDUAL FRAMES -----------')

% Obtain center from the smoothed segmentations (we use the segmentation
% with larger volume (i_meshing)
% Notice that the i_meshing refers to the image indices starting from 0, so
% here corresponds to frame i_meshing+1 in the matlab ordering structure
iLALAA = find(data.img(:,:,:,i_meshing+1) == tags.LA |...
              data.img(:,:,:,i_meshing+1) == tags.LAA);
center = [(max(data.X(iLALAA))+min(data.X(iLALAA)))/2 ...
          (max(data.Y(iLALAA))+min(data.Y(iLALAA)))/2 ...
          (max(data.Z(iLALAA))+min(data.Z(iLALAA)))/2 ];
center = double(center);

if exist('reload_mesh_path')
   if isfile(reload_mesh_path)
      disp(['Loading mesh data from ' reload_mesh_path]);
      load(reload_mesh_path);
   end
else
   % Extract the mesh from the segmentation read at seg_path using only frame
   % i_meshing
   problem_frames =[];
   errors = [];
   tic
   for i = 1:nt_seg
       % Modify faces_cyl_thresh for specific frames and edges (PVs or MV)
       if exist('i_mod','var')
          ii = i_mod == i;
          if sum(ii)>=1
             faces_cyl_thresh = faces_cyl_thresh0.*ones(1,nvein+1);
             faces_cyl_thresh(ied_mod(ii)) = faces_cyl_thresh_mod(ii);
          end
       else
          faces_cyl_thresh = faces_cyl_thresh0;
       end
       disp([' ################## Meshing Frame ' num2str(i) ' ##################']);
       try
          if exist('faces_cyl_thresh','var')
             temp_mesh(i) = mask2mesh(data,i,mesh_keepratio,smooth_iter,MaxArea,...
                                      tags,faces_cyl_thresh,folder_debug);
          else
             temp_mesh(i) = mask2mesh(data,i,mesh_keepratio,smooth_iter,MaxArea,...
                                      tags,0.5,folder_debug);
          end
       catch e
           warning([' ################## PROBLEM WITH Frame ' num2str(i) ' ##################']);
           problem_frames = [problem_frames i];
           errors = [errors e];
       end
   end
   if ~isempty(problem_frames)
       disp('CANNOT CONTINUE! Please fix frames:');
       problem_frames
       for i = 1:length(errors)
           disp(errors(i).identifier);
           disp(errors(i).message);
       end
       return
   end
   time_mesh = toc; 

   % Save data for restart
   vars_to_save = {'iLALAA','center','temp_mesh'};
   save(fullfile(folder_debug,'temp_mesh.mat'),vars_to_save{:},'-v7.3');
end

%% - - - - - - - - REGISTRATION - - - - - - - -

disp(' ----------- REGISTRATION -----------')

if exist('reload_registration_path')
   if isfile(reload_registration_path)
      disp(['Loading registration data from ' reload_registration_path]);
      load(reload_registration_path);
   end
else
   timer_cpd = tic;
   % Displace using the CPD alogrithm
   % Notice that the i_meshing refers to the image indices starting from 0, so
   % here corresponds to frame i_meshing+1 in the matlab ordering structure
   mesh = Mesh_MoveWithCPD(temp_mesh,i_meshing+1,nt_seg,folder_debug);
   toc(timer_cpd);
   mesh = Mesh_AnalyzePM(mesh);

   %% DEBUG
   for it = 1:size(mesh.v_mov,3)
       WriteToVTU(fullfile(folder_debug,...
                           sprintf('debug_mesh_AfterRegistration.%02i.vtu',it)),...
                  mesh.v_mov(:,:,it),mesh.f,'areas',mesh.areas_mov(:,it),...
                  'quality',mesh.qual_mov(:,it),'face_i',mesh.face_type,...
                  'tri_i',1:size(mesh.f,1),'normal_v',mesh.normals_mov(:,1,it),...
                  mesh.normals_mov(:,2,it),mesh.normals_mov(:,3,it));
   end
   %%

   % Save data for restart
   save(fullfile(folder_debug,'temp_registration.mat'),'mesh','-v7.3');
end

%clear temp_mesh;

% Create ib_mesh using triangular faces centers as IB points
[mesh,pp] = Mesh_SmoothSpline(mesh,dt_spln,'seg_time',seg_time,...
                              'ifcubicspln',true);
mesh = Mesh_AnalyzePM(mesh);

%% DEBUG
for it = 1:size(mesh.v_mov,3)
    WriteToVTU(fullfile(folder_debug,...
                        sprintf('debug_mesh_AfterSmoothSpline.%02i.vtu',it)),...
               mesh.v_mov(:,:,it),mesh.f,'areas',mesh.areas_mov(:,it),...
               'quality',mesh.qual_mov(:,it),'face_i',mesh.face_type,...
               'tri_i',1:size(mesh.f,1),'normal_v',mesh.normals_mov(:,1,it),...
               mesh.normals_mov(:,2,it),mesh.normals_mov(:,3,it));
end
%%

% Analyze all the timeframes (needed for next processing steps)
if exist('reload_corks_path')
   if isfile(reload_corks_path)
      disp(['Loading corks data from ' reload_corks_path]);
      load(reload_corks_path);
   end
else
   % Now we have the whole triangular mesh with wall, inlets
   % Let's make the corks
   %mesh = Mesh_AddCorks(mesh,corklen,d_plane);
   % This new routine also moves the corks to avoid them piercing the LA
   mesh = Mesh_AddAndMoveCorks(mesh,corklen,d_plane,mesh.nt);
   mesh = Mesh_AnalyzePM(mesh);

   % Save data for restart
   vars_to_save = {'pp','mesh'};
   save(fullfile(folder_debug,'temp_corks.mat'),vars_to_save{:},'-v7.3');
end

% Center so that LA centered in origin and go from mm to cm
mesh.center = center; % Using centered segmentation
mesh = Mesh_CenterRescale(mesh,0.1);

% Save dt between frames in structure mesh
mesh.dt = T/mesh.nt;

% Smooth triangles' vertices positions in time using Fourier
mesh = Mesh_SmoothFourier2(mesh,nmodes);
mesh = Mesh_AnalyzePM(mesh);

%% DEBUG
for it = 1:size(mesh.v_mov,3)
    WriteToVTU(fullfile(folder_debug,...
                        sprintf('debug_mesh_AfterFourier.%02i.vtu',it)),...
               mesh.v_mov(:,:,it),mesh.f,'areas',mesh.areas_mov(:,it),...
               'quality',mesh.qual_mov(:,it),'face_i',mesh.face_type,...
               'tri_i',1:size(mesh.f,1),'normal_v',mesh.normals_mov(:,1,it),...
               mesh.normals_mov(:,2,it),mesh.normals_mov(:,3,it));
end
%%

disp('Creating Lagrangian mesh of IB points...');
[mesh,ib_mesh] = tri2ibpts2(mesh);

% Adds to the ib_mesh the mesh vertices on the edge of each of the PV.
% Normals and areas are computed approximately since these points don't
% have a triangular faces associated
ib_mesh = IB_AddWallsFromEdges2(ib_mesh,mesh,corklen,0.1*d_plane/2);

%% DEBUG
for it = 1:size(mesh.v_mov,3)
    mesh.qual_mov(:,it) = meshquality(mesh.v_mov(:,:,it),mesh.f);
    WriteToVTU([folder_out 'LA_TriMovingMesh_' num2str(length(mesh.f(:,1))) ...
        'tri_' num2str(it,'%02i') '.vtu'],...
        mesh.v_mov(:,:,it),mesh.f,...
        'areas',mesh.areas_mov(:,it),...
        'quality',mesh.qual_mov(:,it),...
        'face_i',mesh.face_type,....
        'tri_i',1:size(mesh.f,1),...
        'normal_v',mesh.normals_mov(:,1,it),mesh.normals_mov(:,2,it),...
        mesh.normals_mov(:,3,it));

    WriteToVTU([folder_out 'LA_IBMovingMesh_' num2str(length(ib_mesh.type)) ...
        'pts_' num2str(it,'%02i') '.vtu'],...
        ib_mesh.pts(:,:,it),1:1:size(ib_mesh.pts(:,:,it),1),...
        'cell_type',ones([1 size(ib_mesh.pts(:,:,it),1)]),...
        1:1:size(ib_mesh.pts(:,:,it),1),...
        'type_i',ib_mesh.type,'area',ib_mesh.areas(:,it),...
        'vel_v',ib_mesh.vel(:,1,it),ib_mesh.vel(:,2,it),ib_mesh.vel(:,3,it),...
        'normals_v',ib_mesh.normals(:,1,it),ib_mesh.normals(:,2,it),...        
        ib_mesh.normals(:,3,it));    
end
%%

% Check Eulerian bounds set by user
[x0 y0 z0 xf yf zf nx ny nz] = CheckEulerianBounds(x0,y0,z0,xf,yf,zf,...
                                                   nx,ny,nz,mesh,ib_mesh,dh);

% Save data to setup automatic generation of TUCAN tags and launchers
keys = {'caseid','nmodes','resnm','nvein','res','Lc',...
       'x0','xf','y0','yf','z0','zf','dt_spln'};
values = {caseid, nmodes, resnm, nvein, res, Lc,x0,xf,y0,yf,z0,zf,dt_spln};
create_input_file_for_tag_and_launcher_automatic_generation(folder_out,keys,values);

% Create EULERIAN mesh (uniform in 3D):
[xux,yux,zux,...
 xuy,yuy,zuy,...
 xuz,yuz,zuz,...
 xp, yp, zp  ] = mesheu(x0,xf,y0,yf,z0,zf,nx,ny,nz,...
                        'perx',0,'pery',0,'perz',0,...
                        'ghost',true,'write',true);
%
system(['mv mesheu.h5mesh '...
        folder_out sprintf('/%s_%ix%ix%i_eu_v3.h5mesh',casename,nx,ny,nz)]);

% Sort IB points according to type
[type_order,index_order] = sort(ib_mesh.type);
%
pts_order       = ib_mesh.pts(index_order,:,:);
area_order      = ib_mesh.areas(index_order,:);
normal_order    = ib_mesh.normals(index_order,:,:);
vecprod_order   = ib_mesh.vecprod(index_order,:,:);
vel_order       = ib_mesh.vel(index_order,:,:);
% Sort also the fourier coefficient
pts_hat_order     = ib_mesh.pts_hat(index_order,:,:);
area_hat_order    = ib_mesh.areas_hat(index_order,:,:);
normal_hat_order  = ib_mesh.normals_hat(index_order,:,:);
vecprod_hat_order = ib_mesh.vecprod_hat(index_order,:,:);
vel_hat_order     = ib_mesh.vel_hat(index_order,:,:);

% Find the starting and ending indices of each type set
indices_end = find(diff(type_order));
i_start = [1; indices_end+1];
i_end = [indices_end; length(type_order)];
%
nbody = length(i_start(:));
xyz0  = real(squeeze(pts_hat_order(:,:,1)));
vol0  = real(squeeze(area_hat_order(:,1)))'*dh;
normvec0  = real(squeeze(vecprod_hat_order(:,:,1)));

% define more variables needed to create lagrangian geometry with body2file:
ndim=3;
ntensor = 0.5*ndim*(ndim+1);
spi = zeros(ntensor,nbody);
rhoratio=1;
bodyvol=1;
ib = i_start;
ie = i_end;
%
xcontrol = zeros(nbody,3);
for ibody = 1:nbody
 ii = ib(ibody):ie(ibody);
 xcontrol(ibody,:) = mean(xyz0(ii,:));
end

% Create LAGRANGIAN mesh:
fnm = sprintf('%s_%s_lg_v3.dat',casename,resnm);
geometry.body2file(xyz0,vol0,xcontrol,ib,ie,spi,rhoratio,bodyvol,fnm,...
                   'normvec',normvec0,'path',folder_out);

% Fourier modes
fnm = sprintf('%s_%s_fou_v3.dat',casename,resnm);
pts = pts_hat_order(:,:,2:nmodes1+1);
vecprod = vecprod_hat_order(:,:,2:nmodes2+1);
geometry.fou2file(pts,vecprod,nmodes1,nmodes2,fnm,'path',folder_out);

if ifdoublewall
   % Create files (readable by TUCAN) with LAGRANGIAN mesh and fourier modes
   % including a double wall
   %
   % Save values of LAGRANGIAN mesh in a struct (lg)
   lg.prec = 8; lg.ndim = size(xyz0,2); lg.nreal = size(xyz0,1); lg.xyz = xyz0;
   lg.vol = vol0'; lg.normvec = normvec0; lg.xyzc = xcontrol; lg.ib = ib;
   lg.ie = ie; lg.spi = spi'; lg.rhoratio = rhoratio; lg.bodyvol = bodyvol;
   % Save values of fourier modes in a struct (fou)
   fou.ndim = size(xyz0,2); fou.nreal = size(xyz0,1); fou.nm1 = nmodes1;
   fou.nm2 = nmodes2; fou.xyz_r = 2*real(squeeze(pts));
   fou.xyz_i = -2*imag(squeeze(pts)); fou.vecprod_r = 2*real(squeeze(vecprod));
   fou.vecprod_i = -2*imag(squeeze(vecprod));
   % Generate new lg and fou structs
   nveins = mesh.n_edges-1;
   [lg_dw,fou_dw] = geometry.gen_doubleW(lg,fou,nveins); clear lg fou;
   % Save LAGRANGIAN mesh file including a double wall
   lg_fnm = sprintf('%s_%s_lg_dw.dat',casename,resnm);
   geometry.body2fileS(lg_dw,lg_fnm,'path',folder_out);
   % Save fourier modes file including a double wall
   fou_fnm = sprintf('%s_%s_fou_dw.dat',casename,resnm);
   geometry.fou2fileS(fou_dw,fou_fnm,'path',folder_out);
end

% Calculate variables required to create LA mask in postprocessing
caras = mesh.f(mesh.face_type<=mesh.n_edges,:);
coeffs = mesh.v_hat;
norm_out = mesh.norm_out;
% Create .mat files with data to create LA mask in postprocessing
filename = fullfile(folder_out,sprintf('data_for_mask_%ix%ix%i',nx,ny,nz));
save(filename,'caras','coeffs','nmodes1','norm_out');

% Save LAA mask
LAAarea = LAAmask_generator(folder_out,tags.LAA,data,mesh,nmodes,nmodes1,...
                            i_meshing,T,nx,ny,nz,'seg_time',seg_time);

% Create Qmifile.dat for TUCAN
seg2QmitFile_fun(seg_path,seg_ext,tags.LV,T,nmodes,folder_out,caseid,...
                 'seg_time',seg_time);

% Move data saved during registration and debugging in folder_out
disp(sprintf('Moving data from:\n %s\nto:\n %s',folder_debug,folder_out));
mkdir([folder_out 'cpd_registration/']);
movefile([folder_debug 'CPD_dump*.mat'],[folder_out 'cpd_registration/']);
mkdir([folder_out 'single_meshes_smooth_edges/']);
movefile([folder_debug 'debug_mesh.*.vtu'],...
         [folder_out 'single_meshes_smooth_edges/']);
mkdir([folder_out 'debug/']);
mkdir([folder_out 'debug/registered_meshes/']);
movefile([folder_debug 'debug_mesh_AfterRegistration.*.vtu'],...
         [folder_out 'debug/registered_meshes/']);
mkdir([folder_out 'debug/smooth_splined_meshes_with_corks/']);
movefile([folder_debug 'debug_mesh_AfterSmoothSpline.*.vtu'],...
         [folder_out 'debug/smooth_splined_meshes_with_corks/']);
mkdir([folder_out 'debug/smoothed_splined_meshes_with_corks_and_smooth_edges/']);
movefile([folder_debug 'debug_mesh_AfterFourier.*.vtu'],...
         [folder_out 'debug/smooth_splined_meshes_with_corks/']);
mkdir([folder_out 'debug/SmoothFourier2/']);
movefile([folder_debug 'LAmesh_full_prefix*.vtu'],...
         [folder_out 'debug/SmoothFourier2/']);
movefile([folder_debug 'data_levset.mat'],[folder_out 'debug/']);
movefile([folder_debug 'temp_*.mat'],[folder_out 'debug/']);

% Save whole workspace
save([folder_out 'full_preproc.mat'],'-v7.3');

% Check mesh element size
fs1 = 14; fs2 = 18;
figelesiz = figure('Visible','off'); 
hist(sqrt(mesh.areas_mov(:,i_meshing+1)));
yl = ylim; hold on; plot([dh dh],yl,'-r','LineWidth',2);
set(gca,'FontSize',fs1,'TickLabelInterpreter','Latex');
xlabel('triangle characteristic length (cm)','FontSize',fs2,'Interpreter','Latex');
ylabel('Number of triangles','FontSize',fs2,'Interpreter','Latex');
box on; grid on;
figelesiznm = fullfile(folder_out,'larger_mesh_elements_size.png');
print(figelesiz,figelesiznm,'-dpng','-r300');
toc(timer_start)

exit;
