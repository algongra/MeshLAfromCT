function data = CT_GetSmoothSegWithLevSet(seg_path, nt_seg, seg_ext, tags,...
    smoothing_stencil, filter_type)

if nargin<6
    filter_type = 'avg';
end

% Compute level set to smooth in space and time
[lev_set, data] = CT_ComputeLevSet(seg_path, nt_seg, seg_ext, [tags.LA tags.LAA]);
data.img = zeros([data.ny data.nx data.nz nt_seg],'int8');
data.img = CT_AddLevSetWithTag(lev_set, data.img, tags.LA, smoothing_stencil, filter_type);
% WriteToVTR('segLA.vtr',1:size(data.img,1),1:size(data.img,2),1:size(data.img,3),'seg',data.img(:,:,:,1));

lev_set = CT_ComputeLevSet(seg_path, nt_seg, seg_ext, tags.LAA);
data.img = CT_AddLevSetWithTag(lev_set, data.img, tags.LAA-tags.LA, smoothing_stencil, filter_type);
% WriteToVTR('segLAA.vtr',1:size(data.img,1),1:size(data.img,2),1:size(data.img,3),'seg',data.img(:,:,:,1));

lev_set = CT_ComputeLevSet(seg_path, nt_seg, seg_ext, tags.LV);
data.img = CT_AddLevSetWithTag(lev_set, data.img, tags.LV, smoothing_stencil, filter_type);
% WriteToVTR('segLV.vtr',1:size(data.img,1),1:size(data.img,2),1:size(data.img,3),'seg',data.img(:,:,:,1));

for tag = tags.PVs
    lev_set = CT_ComputeLevSet(seg_path, nt_seg, seg_ext, tag);
    % WriteToVTR(['lsPV' num2str(tag) '.vtr'],...
    %    1:size(lev_set,1),1:size(lev_set,2),1:size(lev_set,3),'ls',lev_set(:,:,:,1)); 
    data.img = CT_AddLevSetWithTag(lev_set, data.img, tag, smoothing_stencil, filter_type);
    % WriteToVTR(['segPV' num2str(tag) '.vtr'],...
    %    1:size(data.img,1),1:size(data.img,2),1:size(data.img,3),'seg',data.img(:,:,:,1));
end

return
