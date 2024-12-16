function mesh_out = CT_MeshAndResample(data,i,mesh_keepratio,smooth_iter,hstencil)

if nargin<4
    smooth_iter = 2;
end

if nargin<5
    hstencil = [3 3 3];
end

% Added for legacy:
if length(hstencil) == 1
    hstencil = hstencil*[1 1 1];
end

disp(sprintf('\nObtaining initial mesh...'));

if any(hstencil > 1)
    h3 = ones(hstencil);
    h3 = h3./ sum(h3(:));
    mask = imfilter(data.mask(:,:,:,i),h3);
    mask = imfill(mask>0.5,'holes');
else
    mask = single(data.mask(:,:,:,i));
end

[f,v] = isosurface(data.X,data.Y,data.Z,mask);

%%%%%% DEBUG Checks normals now
% normals = surfacenorm(v,f);
% WriteToVTU('debug_isosurf.vtu',v,f,'normal_v',normals(:,1),...
%     normals(:,2),normals(:,3));
%%%%%%

if mesh_keepratio<1
    disp(sprintf('\nResampling Mesh...'));
    % Resample faces
    [v,f]=meshresample(v,f,mesh_keepratio);
end

%%%%%% DEBUG Checks normals now
% normals = surfacenorm(v,f);
% WriteToVTU('debug_resample.vtu',v,f,'normal_v',normals(:,1),...
%     normals(:,2),normals(:,3));
%%%%%%

if mesh_keepratio >= 0.1 && smooth_iter > 0
    % Smooth the surface only if we keep many elements, otherwise it should
    % be already smooth enough
    disp(sprintf('\nSmoothing Mesh...'));
    % f = unique(sort(f,2),'rows'); % THIS MESSES UP THE ORDER AND THE
    % NORMALS: DO NOT USE!!!
    [conn,connnum,count] = meshconn(f,length(v));
    % v = smoothsurf(v,[],conn,smooth_iter,1,'lowpass');
    v = smoothsurf(v,[],conn,smooth_iter,1);
end
% %%%%%% DEBUG Checks normals now
% normals = surfacenorm(v,f);
% WriteToVTU('debug_smooth.vtu',v,f,'normal_v',normals(:,1),...
%     normals(:,2),normals(:,3));
% %%%%%%

% Use structure to output
mesh_out.v = v;
mesh_out.f = f;
mesh_out.keep_ratio = mesh_keepratio;
return
