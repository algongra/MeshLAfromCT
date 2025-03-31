function WriteSmoothedNifti(seg_path,seg_ext,data,folder_debug);

%%%%% paths for additional routines
% Load NIfTI images 
addpath('../NIfTI_load');

% Read image header and aditional info to copy it in the new image
seg_fullfile = fullfile(seg_path,sprintf('0%s',seg_ext));
nii_data = load_nii(seg_fullfile);

% Copy header and additional info to new image
nii.hdr = nii_data.hdr;
nii.original = nii_data.original;
% Create new info (original struct is not modified to keep information)
%
% Modify hdr.dime.struct
nii.hdr.dime.dim(1:5) = [ndims(data.img) size(data.img,1) size(data.img,2) ...
                         size(data.img,3) size(data.img,4)];
nii.hdr.dime.pixdim(2:4) = [data.dx data.dy data.dz];
% Modify hist.dime.struct
nii.hdr.hist.srow_x(4) = data.dx;
nii.hdr.hist.srow_y(4) = data.dy;
nii.hdr.hist.srow_z(4) = data.dz;
% Modify machine field
nii.machine = 'Smoothed_Segmentation_By_MeshLAfromCT';
% Copy image from smoothed data
nii.img = data.img;
% Define filename of smoothed segmentation
filename = fullfile(folder_debug,'smoothed_segmentation.nii.gz');

% Save Smoothed segmentation
save_nii(nii, filename);

return

end
