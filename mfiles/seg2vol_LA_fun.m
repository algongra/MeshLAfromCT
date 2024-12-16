function i_meshing = seg2vol_LA_fun(seg_path,seg_ext,seg_tags,seg_nt,T);
% @author A.Gonzalo
%
% @brief Funcion to print LA volume info and to obtain i_meshing from
%        segmentations files
%
% @date 03-30-2023 by A.Gonzalo \n
%                  Created and documented
% @details
%
% MANDATORY ARGUMENTS
% -------------------
%  - seg_path: path where segmentation files ([0-20].nii or [0-20].nii.gz) are
%              stored. [char]
%  - seg_ext:  extension of the segmentation files (.nii, or .nii.gz). [char]
%  - seg_tags: tags or labels identifying the LA body, and the LAA in the
%              segmentation files (tags MUST be equal in all segmentation files)
%              [single, double, int8,...]. Size: 1x2
%  - T:        cardiac cycle duration in seconds (s) [single or double].
%              Size: 1x1
%
% OUTPUT
% ------
%  - i_meshing: index of segmentation file where LA volume is maximum


% Load library to read NIfTI images
%addpath('ct_processing_la-20240109T001842Z-001\ct_processing_la\NIfTI_load');
addpath('../NIfTI_load');

% Initialize array
vol_vs_time_vox = zeros([1 seg_nt]);
fprintf('\nReading segmentation files........');
for i = 1:seg_nt
    fprintf('\b\b\b\b\b%02i/%02i',i,seg_nt);
    
    % Load file
    seg_fullfile = fullfile(seg_path,[num2str(i-1) seg_ext]);
    nii_data = load_nii(seg_fullfile);
    
    % Volume by summing pixels of the different tags
    for itag = 1:length(seg_tags)
        vol_vs_time_vox(i) = vol_vs_time_vox(i) +...
                             sum(nii_data.img(:)==seg_tags(itag));
    end
end

% Get voxel size in [cm]
dx = nii_data.hdr.dime.pixdim(2)/10;
dy = nii_data.hdr.dime.pixdim(3)/10;
dz = nii_data.hdr.dime.pixdim(4)/10;

clear nii_data;

% Volume in [cm^3]
vol_vs_time = vol_vs_time_vox.*dx.*dy.*dz;

% time span
time = linspace(0,T,length(vol_vs_time)+1); time = time(1:end-1);

% Plot volume vs time
fs1 = 14; fs2 = 18;
plot(time,vol_vs_time,'-','LineWidth',2,'Color','b'); hold on;
set(gca,'FontSize',fs1,'TickLabelInterpreter','Latex');
xlabel('time (s)','FontSize',fs2,'Interpreter','Latex');
ylabel('$\mathrm{Volume}_{\mathrm{LA}}$ $(\mathrm{cm}^3)$','FontSize',fs2,...
       'Interpreter','Latex');
box on; grid on; xlim([0 T]);

% Obtain and print i_meshing, maximum and minimum LA volume, and LA Emtying
% Fraction (EF) based on the segmentations (not interpolations performed)
[LAvolmax,i_meshing] = max(vol_vs_time);
i_meshing = i_meshing - 1; % i_meshing is an image segmetation index
                           % (segmentation indices start with 0, not with 1)
LAvolmin = min(vol_vs_time);
%
disp(sprintf('\n\ni_meshing = %i\n',i_meshing));
%
disp(sprintf('Based on the segmentations (not interpolation performed):'));
disp(sprintf(' - Minimum LA volume = %5.2f',LAvolmin));
disp(sprintf(' - Maximum LA volume = %5.2f',LAvolmax));
disp(sprintf(' - LA emptying fraction (EF) = %4.2f\n',1-LAvolmin/LAvolmax));

return;

end

