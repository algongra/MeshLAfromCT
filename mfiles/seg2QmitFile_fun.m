function seg2QmitFile_fun(seg_path,seg_ext,seg_tags,T,nmodes,out_path,caseid,varargin);
% @author A.Gonzalo
%
% @brief Funcion to obtain the flow through the Mitral valve and write a
%        TUCAN readable file with this information (required to run simulations)
%
% @date 08-11-2023 by A.Gonzalo \n
%                  Created and documented
% @details
%
% MANDATORY ARGUMENTS
% -------------------
%  - seg_path: path where segmentation files ([0-20].nii or [0-20].nii.gz) are
%              stored. [char]
%  - seg_ext:  extension of the segmentation files (.nii, or .nii.gz). [char]
%  - seg_tags: tag(s) or label(s) identifying the LV in the segmentation files
%              (tag(s) number MUST be the same in all segmentation files)
%              [integer: single, double, int8,...]. Size: 1x1
%  - T:        cardiac cycle period in [s]. [single, double]. Size: 1x1
%  - nmodes:   number of Fourier modes to reconstruct the Qmit vs time signal.
%              [integer: single, double, int8,...]. Size: 1x1
%  - out_path: path where function stores output file. [char]
%  - caseid:   ID of subject analyzed. [char]
%
% OPTIONAL ARGUMENTS
% ------------------
%  - seg_time:   time instants at which segmentation files (CT-scan DICOMs)
%                were obtained along the cardiac cycle in [s]. [double]. 
%                Size: 1xnt, where nt is the number of segmentations used
%                ACHTUNG: seg_time MUST NOT have time instants from two
%                         different cardiac cycles (e.g., 0 and T MUST NOT be
%                         included, only one of those two values MUST be
%                         included). Default: []
%    dt_sim:     dt of nominal spatial resolution simulation in TUCAN.
%                [single, double]. Size: 1x1. Default: 5e-5
%  - Qmitflnm:   beginning of Qmit OUTPUT file name (readable by TUCAN). [char].
%                Default: 'Qmit_vs_t_one_cycle'
%  - dt:         dt to interpolate dvdt_s signal obtained from Fourier into Qmit
%                written in Qmitflnm. [single, double]. Size: 1x1.
%                Default: 0.001
%  - erasewaves: array identifying the waves to be deleted.
%                waves are plotted in yellow in $Q_mit$ vs $time/T figure when
%                this function is called.
%                waves id increases from left to right starting by 1, until n.
%                [integer: single, double, int8,...]. Size: 1x1. Default: []
%    csapsparamsdefault: Use default weights parameters of csaps function
%                        [logical] Size: 1x1. Default: False
%
% OUTPUT
% ------
%  - File ${caseid}_vol_vs_time.mat created in folder ${out_path}.
%  - File ${Qmitflnm}_${caseid}.dat created in folder ${out_path}.
%    This file is an input required by TUCAN to calculate the boundary
%    conditions at the pulmonary veins [PVs] (i.e., the velocity to be
%    prescribed at the PVs inlets).
%
% EXAMPLES
% --------
%  @code
%  seg_path = '/Users/alex/Projects/Hearts/Bolus/CT_scans_data/26393/seg-nii_LV';
%  seg_ext = '.nii.gz';
%  seg_tags = 1;
%  T = 1;
%  nmodes = 5;
%  out_path = '/Users/alex/Projects/Hearts/ct_processing_la/26393_5modes_lowres_0.2_20230518_1337';
%  caseid = '26393';
%  seg2QmitFile_fun(seg_path,seg_ext,seg_tags,T,nmodes,out_path,caseid);
%  @endcode
%
%  @code
%  seg_path = '/Users/alex/Projects/Hearts/Bolus/CT_scans_data/26393/seg-nii_LV';
%  seg_ext = '.nii.gz';
%  seg_tags = 1;
%  T = 1;
%  nmodes = 10;
%  out_path = '/Users/alex/Projects/Hearts/ct_processing_la/26393_5modes_lowres_0.2_20230518_1337';
%  caseid = '26393';
%  seg2QmitFile_fun(seg_path,seg_ext,seg_tags,T,nmodes,out_path,caseid);
%  @endcode
%
%  @code
%  seg_path = '/users/alex/projects/hearts/bolus/ct_scans_data/28844/seg-nii_LV';
%  seg_ext = '.nii.gz';
%  seg_tags = 1;
%  T = 1;
%  nmodes = 5;
%  out_path = '/users/alex/projects/hearts/ct_processing_la/28844_5modes_nomres_0.4_20230517_1649';
%  caseid = '28844';
%  seg_time = [0 0.1 0.2 0.3 0.35 0.4 0.5 0.6 0.7 0.8 0.9];
%  seg2QmitFile_fun(seg_path,seg_ext,seg_tags,T,nmodes,out_path,caseid,...
%                   'seg_time',seg_time);
%  @endcode
%
%  @code
%  seg_path = '/users/alex/projects/hearts/bolus/ct_scans_data/28844/seg-nii_lv';
%  seg_ext = '.nii.gz';
%  seg_tags = 1;
%  T = 1;
%  nmodes = 5;
%  out_path = '/users/alex/projects/hearts/ct_processing_la/28844_5modes_nomres_0.4_20230517_1649';
%  caseid = '28844';
%  qmitflnm = 'qmit_vs_t_deleting_first_wave';
%  waveid = [1];
%  seg2QmitFile_fun(seg_path,seg_ext,seg_tags,T,nmodes,out_path,caseid,...
%                   'qmitflnm',qmitflnm,'erasewaves',waveid);
%  @endcode

% Load library to read NIfTI images
addpath('../NIfTI_load');
% Load IBcode Matlab library (to use misc.assigndefaults function)
addpath('~/IBcode/lib/mfiles');


% defaults
seg_time = [];
dt_sim = 5e-5;
Qmitflnm = 'Qmit_vs_t_one_cycle';
dt = 0.001;
erasewaves = [];
csapsparamsdefault = false;
misc.assigndefaults(varargin{:});


% Check number of segmentations saved in seg_path directory
filelist = dir(fullfile(seg_path,['*' seg_ext]));
% Get number of segmentations available (using filelist)
nt_seg = length(filelist);
% Sanity checks
if nt_seg == 0
   error(sprintf('No segmentation files found in\n%s',seg_path));
end
%
if ~isempty(seg_time)
   % Get number of segmentation (using seg_time)
   seg_nt = length(seg_time);
   % Sanity check
   if numel(seg_nt)~=numel(nt_seg)
      error(sprintf([' seg_time array MUST have the same number of elements '...
                     'as number of segmentation files stored in\n%s'],seg_path));
   end
   % Create time span by normalizing seg_time with the cardiac cycle period and
   % offseting it seg_time(1) to start the cycle at 0
   time = (seg_time - seg_time(1))./T;
else
   % Create time span with constant dt
   % (if seg_time not defined by user as varargin)
   time = linspace(0,T,nt_seg+1); time = time(1:end-1);
end

% Initialize array to save LV volume over time in [voxels]
vol_vs_time_vox = zeros([1 nt_seg]);
fprintf('\nReading segmentation files........');
for i = 1:nt_seg
    fprintf('\b\b\b\b\b%02i/%02i',i,nt_seg);
    
    % Load segmentation file
    seg_fullfile = fullfile(seg_path,[num2str(i-1) seg_ext]);
    if isfile(seg_fullfile)
       nii_data = load_nii(seg_fullfile);
    else
       error(sprintf(' File:\n %s\n does not exist',seg_fullfile));
    end
    
    % Calculate LV volume from this segmentation file by summing pixels of the
    % different tags
    for itag = seg_tags
        vol_vs_time_vox(i) = vol_vs_time_vox(i) +...
                             sum(nii_data.img(:)==seg_tags(itag));
    end
end

% Get voxel size in [cm]
dx = nii_data.hdr.dime.pixdim(2)/10;
dy = nii_data.hdr.dime.pixdim(3)/10;
dz = nii_data.hdr.dime.pixdim(4)/10;

clear nii_data;

% Obtain LV volume over volume in [cm^3]
vol_vs_time = vol_vs_time_vox.*dx.*dy.*dz;

% Create time span required in nominal spatial resolution simulation in TUCAN
t_sim = [0:dt_sim:1];

% Check if dt is constant among the normalized time instants
isdtcte = abs(sum(diff(time)-time(2)))<1e-10;

% Obtain dvdt signal (flow through the MV) using Fourier
%
% If enough segmentations are available for the nmodes Fourier selected and if
% dt is constant
if nmodes<ceil(nt_seg/2) & isdtcte
   time_s = time;
   % Time/frequency analysis
   dt_aux = mean(diff(time_s));
   k   = 2*pi*(-1/(2*dt_aux):(1/(dt_aux*nt_seg)):(1/(2*dt_aux)));
   kt  = ifftshift(k(1:end-1));

   % Fourier transform of LV volume signal obtained from segmentations
   vol_vs_time_hat = fft(vol_vs_time);

   % Derive LV volume in Fourier space to obtain dvdt
   dvdt_hat = 1i.*kt.*vol_vs_time_hat;
   dvdt = real(ifft(dvdt_hat));

   % Drop some frequencies in _s arrays
   % LV volume
   vol_vs_time_hat_s = vol_vs_time_hat;
   vol_vs_time_hat_s(round(nmodes)+1:end-(round(nmodes)-1)) = 0;
   vol_vs_time_s = real(ifft(vol_vs_time_hat_s));
   % dvdt (flow through the MV)
   dvdt_hat_s = 1i.*kt.*vol_vs_time_hat_s;
   dvdt_s = real(ifft(dvdt_hat_s));

   % Plot
   fs1 = 14; fs2 = 18;
   % LV volume obtained with Fourier vs time
   figure();
   plot(time,vol_vs_time_s,'--','LineWidth',2,'Color',[0 0.5 0]); hold on;
   scatter(time,vol_vs_time,100,'ob','MarkerEdgeColor','none',...
           'MarkerFaceColor','b','MarkerFaceAlpha',0.5);
   set(gca,'FontSize',fs1,'TickLabelInterpreter','Latex');
   xlabel('$time/T$','FontSize',fs2,'Interpreter','Latex');
   ylabel('$\mathrm{Volume}_{\mathrm{LV}}$ $(\mathrm{cm}^3)$','FontSize',fs2,...
          'Interpreter','Latex');
   box on; grid on; xlim([0 1]);
   % dvdt (flow through MV) obtained with Fourier vs time
   figure();
   plot(time,dvdt,'-.r','LineWidth',2); hold on;
   plot(time,dvdt_s,'--','LineWidth',2,'Color',[0 0.5 0]);
   set(gca,'FontSize',fs1,'TickLabelInterpreter','Latex');
   xlabel('$time/T$','FontSize',fs2,'Interpreter','Latex');
   ylabel('dVoldt $(\mathrm{cm}^3/s)$','FontSize',fs2,'Interpreter','Latex');
   box on; grid on; xlim([0 1]);

% If NOT enough segmentations are available for the nmodes Fourier selected
% and/or if dt is NOT constant
else
   % Get cubic smoothing spline to reconstruct LV volume from vol_vs_time
   %
   % 1) Repeat some points of the vol_vs_time at the beginning at the end to
   %    avoid problem with the periodicity
   n_pts_rep = nt_seg-1;
   vol = [vol_vs_time(end-n_pts_rep:end) vol_vs_time vol_vs_time(1:n_pts_rep)];
   tim = [time(end-n_pts_rep:end) T+time 2*T+time(1:n_pts_rep)];
   if ~csapsparamsdefault
      % 2) Define p (weight between error and roughness measures) and $\lambda$
      %    (weights in the roughness measure)
      %     - p is the first element of 1D array p
      %     - $\lambda$ are 2:end elements of 1D array p
      %       $\lambda = \left(\frac{\Delta Vol_{LV}}{\Delta t}\right)$ [1]
      %    2.1) Get $\lambda$ as in [1]
      n = abs( (vol(2:end) - vol(1:end-1))./(tim(2:end) - tim(1:end-1)) );
      %    2.2) Normalized $\lambda$ ensuring elements are non-zero
      n = max(n/max(n),1e-10); % avoid zeros
      %    2.3) Introduce p at the beginning of 1D array p
      p = [-1 n.^(-1)];
      % 3) Define $\omega$ (weights in the error measure)
      %    - $\omega = \left| Vol_{LV} - mean\left( Vol_{LV} \right)\right|$ [2]
      %    3.1) Start off with the trapezoidal error weights
      dt_seg = diff(tim);
      weights = ([dt_seg 0]+[0 dt_seg])/2;
      %    3.2) get $\omega$ as in [2]
      w = max(weights.*abs(vol - mean(vol)), 1e-10); % avoid zeros
      % 4) No sampling time yet (we will use fnval later)
      xx = [];
   else
      % Start off with the trapezoidal error weights
      dt_seg = diff(tim);
      weights = ([dt_seg 0]+[0 dt_seg])/2;
      % csaps chooses the smoothing parameter (since p = -1)
      p = -1;
      % No sampling time yet (we will use fnval later)
      xx = [];
   end
   % 5) Use csaps to obtain cubic smoothing spline
   spln = csaps(tim,vol,p,xx,weights);
   %
   % Keep the splines required to reconstruct one cycle
   pp.coefs = spln.coefs(1+n_pts_rep:1+n_pts_rep+nt_seg,:);
   pp.breaks = spln.breaks(1+n_pts_rep:1+n_pts_rep+nt_seg+1)-1;
   pp.dim = spln.dim;
   pp.form = spln.form;
   pp.order = spln.order;
   pp.pieces = size(pp.coefs,1);
   %
   % Evaluate splines at all cycle times solved in a simulation
   vspln = fnval(pp,t_sim);
   %
   % Print error in periodicity of LV volume obtained with splines
   disp(sprintf(['\nPeriodicity error in LV volume obtained with splines = '...
                 '%g'],abs(vspln(1)-vspln(end))));

   % Get splines of dvdt (flow through the MV)
   ppdt.coefs = [3*pp.coefs(:,1) 2*pp.coefs(:,2) pp.coefs(:,3)];
   ppdt.breaks = spln.breaks(1+n_pts_rep:1+n_pts_rep+nt_seg+1)-1;
   ppdt.dim = spln.dim;
   ppdt.form = spln.form;
   ppdt.order = spln.order-1;
   ppdt.pieces = size(ppdt.coefs,1);
   %
   % Evaluate splines at all cycle times solved in a simulation
   dvdtspln = fnval(ppdt,t_sim);
   %
   % Print error in periodicity of dvdt (flow through the MV) obtained with
   % splines
   disp(sprintf(['Periodicity error in dLVVoldt (flow through the MVs) '...
                 'obtained with splines = %g'],abs(dvdtspln(1)-dvdtspln(end))));

   % Fourier fitting of smoothed and periodic data from splines
   %
   time_s = t_sim(1:end-1);
   % From volume signal obtained with splines
   vol_vs_time_hat = fft(vspln(1:end-1));
   vol_vs_time_hat_s = vol_vs_time_hat;
   vol_vs_time_hat_s(round(nmodes)+1:end-(round(nmodes)-1)) = 0;
   vol_vs_time_s = real(ifft(vol_vs_time_hat_s));
   %
   % From volume derivative signal obtained with splines
   dvdt_hat = fft(dvdtspln(1:end-1));
   dvdt_hat_s = dvdt_hat;
   dvdt_hat_s(round(nmodes)+1:end-(round(nmodes)-1)) = 0;
   dvdt_s = real(ifft(dvdt_hat_s));
   %{
   %
   % Deriving LV volume signal obtained with splines in Fourier space
   %  - This option yields really similar results than the previous one
   %
   % Time/frequency analysis
   nt_sim = length(t_sim)-1;
   k   = 2*pi*(-1/(2*dt_sim):(1/(dt_sim*nt_sim)):(1/(2*dt_sim)));
   kt  = ifftshift(k(1:end-1));
   %
   dvdt_hat = 1i.*kt.*vol_vs_time_hat;
   dvdt_hat_s = dvdt_hat;
   dvdt_hat_s(round(nmodes)+1:end-(round(nmodes)-1)) = 0;
   dvdt_s = real(ifft(dvdt_hat_s));
   %}

   % Plot
   fs1 = 14; fs2 = 18;
   % LV volume obtained with Fourier vs time
   figure();
   plot(t_sim,vspln,'-k','LineWidth',2); hold on;
   plot(time_s,vol_vs_time_s,'--','LineWidth',2,'Color',[0 0.5 0]);
   scatter(time,vol_vs_time,100,'ob','MarkerEdgeColor','none',...
           'MarkerFaceColor','b','MarkerFaceAlpha',0.5);
   set(gca,'FontSize',fs1,'TickLabelInterpreter','Latex');
   xlabel('$time/T$','FontSize',fs2,'Interpreter','Latex');
   ylabel('$\mathrm{Volume}_{\mathrm{LV}}$ $(\mathrm{cm}^3)$','FontSize',fs2,...
          'Interpreter','Latex');
   box on; grid on; xlim([0 1]);
   % dvdt (flow through MV) obtained with Fourier vs time
   figure();
   plot(t_sim,dvdtspln,'-k','LineWidth',2); hold on;
   plot(time_s,dvdt_s,'--','LineWidth',2,'Color',[0 0.5 0]);
   set(gca,'FontSize',fs1,'TickLabelInterpreter','Latex');
   xlabel('$time/T$','FontSize',fs2,'Interpreter','Latex');
   ylabel('dVoldt $(\mathrm{cm}^3/s)$','FontSize',fs2,'Interpreter','Latex');
   box on; grid on; xlim([0 1]);
end

% Save LV volume over time and Qmit over time signals in
% ${caseid}_vol_vs_time.mat file
save(fullfile(out_path,strcat(caseid,'_vol_vs_time.mat')));


% Preparing data to write file ${Qmitflnm}_${caseid}.dat
fig = figure();

% Sanity check
if dt_sim > dt
   error(sprintf(['dt selected is larger than the dt of the simulation\n'...
                  'subroutine Qmitral_vs_time in user_moveib.f requires '...
                  'dt_sim < dt to work']));
end

% Save time and dvdt signals from Fourier fitting of smoothed and periodic data
% from spline in temporal arrays
ttmp = time_s;
qtmp = dvdt_s;

% Closing the cardiac cycle with the first value (it is periodic in time)
yy = [qtmp qtmp(1)];
xx = [ttmp T];

% Interpolate data to new values of time defined with the dt selected
tt = [xx(1):dt:xx(end)];
Qmit = spline(xx,yy,tt);

% Plot interpolated data
set(0,'CurrentFigure',fig); plot(tt,Qmit,'-b','LineWidth',2); hold on;

% Remove negative differential of volume (when LV contracts, MV is closed and
% no flow returns to the LA)
Qmit = max(Qmit,0);

% Plot Qmit once negative values have been removed
set(0,'CurrentFigure',fig); plot(tt,Qmit,'--c','LineWidth',2);

% Calculate mask to ease waves visualization
mask = diff(Qmit)~=0;
mask = [0 mask];
iwave = find(diff(mask)~=0);
maskappend = mask(iwave); %iwave = iwave + maskappend;
Qmitaux = Qmit; ttaux = tt;
icnt = 0;
for ii = 1:length(iwave)
    icnt = icnt + 1;
    iwv = iwave(ii);
    if maskappend(icnt)
       mask = [mask(1:iwv) 0 mask(iwv+1:end)];
       Qmitaux = [Qmitaux(1:iwv) Qmitaux(iwv-1) Qmitaux(iwv+1:end)];
       ttaux = [ttaux(1:iwv) ttaux(iwv) ttaux(iwv+1:end)];
    else
       mask = [mask(1:iwv) 1 mask(iwv+1:end)];
       Qmitaux = [Qmitaux(1:iwv) Qmitaux(iwv+1) Qmitaux(iwv+1:end)];
       ttaux = [ttaux(1:iwv) ttaux(iwv) ttaux(iwv+1:end)];
    end
    iwave = iwave + 1;
end

% Plot waves
set(0,'CurrentFigure',fig);
plot(ttaux,mean(Qmitaux)*mask,'LineWidth',2,'Color',[0.9290 0.6940 0.1250]);

% Remove waves if erasewaves array is defined by user as varargin
if ~isempty(erasewaves)
   waves = bwconncomp(mask);
   % Sanity check
   if any(erasewaves > waves.NumObjects)
      error(sprintf(['erasewaves cannot erase waves with id number lager '...
                     'than the maximum number of waves\nIn this case # of '...
                     'waves = %i'],waves.NumObjects));
   else
      for iw = erasewaves
          Qmit(waves.PixelIdxList{iw}) = 0;
          disp(sprintf('Wave #%i removed',iw));
      end
      disp(sprintf(['Do you still want the output file to begin with %s?'...
                    '\nIf not, remember to also define the varargin Qmitflnm'],...
                    Qmitflnm));
   end
end

% Plot Qmit signal to be saved in ${Qmitflnm}_${caseid}.dat file
set(0,'CurrentFigure',fig); plot(tt,Qmit,'LineWidth',2,'Color',[0 0.5 0]);
set(gca,'FontSize',fs1,'TickLabelInterpreter','Latex');
xlabel('$time/T$','FontSize',fs2,'Interpreter','Latex');
ylabel('$Q_{mit}$ $(\mathrm{cm}^3/s)$','FontSize',fs2,'Interpreter','Latex');
box on; grid on; xlim([0 T]);

% Writing Qmit to file ${Qmitflnm}_${caseid}.dat
nt = length(Qmit(:));
fnm = fullfile(out_path,strcat(Qmitflnm,'_',caseid,'.dat'));
%
fid = fopen(fnm,'w','l');
fwrite(fid,4,'int');
fwrite(fid,nt,'int');
fwrite(fid,4,'int');
%
fwrite(fid,8*nt,'int');
fwrite(fid,tt,'double');
fwrite(fid,8*nt,'int');
%
fwrite(fid,8*nt,'int');
fwrite(fid,Qmit,'double');
fwrite(fid,8*nt,'int');
%
fclose(fid);

return

end

