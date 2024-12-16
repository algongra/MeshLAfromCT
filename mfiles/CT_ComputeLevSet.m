function [lev_set, data] = CT_ComputeLevSet(seg_path, nt_seg, seg_ext, include_tags)
% Extract a level set function from the segmentation, with positive values
% inside of the tagged volume and negative outside.
if nargin < 4
    % Inclusion tags
    include_tags = [4 10];    
end
fprintf(1,'\nComputing level set.');
tic
% Read first frame for initialization
i_seg = 1;
seg_fullfile = fullfile(seg_path,[num2str(i_seg-1) seg_ext]);
% Get volumetric dataset from nifti file
data = CT_nii2vol(seg_fullfile, true);
lev_set = zeros([size(data.img) nt_seg], 'single');
inside_mask = false(size(data.img));
for tag = include_tags
    inside_mask = inside_mask | (data.img==tag);
end
lev_set(:,:,:,i_seg) = single(inside_mask.*bwdist(~inside_mask)-(~inside_mask).*bwdist(inside_mask));

for i_seg = 2:nt_seg
    fprintf(1,'.');
    seg_fullfile = fullfile(seg_path,[num2str(i_seg-1) seg_ext]);
    data = CT_nii2vol(seg_fullfile, true);
    inside_mask = false(size(data.img));
    for tag = include_tags
        inside_mask = inside_mask | (data.img==tag);
    end
    lev_set(:,:,:,i_seg) = single(inside_mask.*bwdist(~inside_mask)-(~inside_mask).*bwdist(inside_mask));
end
data = rmfield(data,'img');
toc
disp(' ');
return
