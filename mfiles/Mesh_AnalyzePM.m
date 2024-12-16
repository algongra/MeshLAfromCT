function mesh = Mesh_AnalyzePM(mesh,it)
% Compute basic mesh properties: areas of face, mesh quality and normals

if ~isfield(mesh,'norm_out')
    mesh = Mesh_CheckNormalsOrientation(mesh);
end

if isfield(mesh,'v_mov') && nargin<2
    % Analyze all meshes
    fprintf('\nAnalyzing all meshes........');    
    nt = size(mesh.v_mov,3);
    if isfield(mesh,'qual_mov'); mesh = rmfield(mesh,'qual_mov'); end;
    if isfield(mesh,'areas_mov'); mesh = rmfield(mesh,'areas_mov'); end;
    if isfield(mesh,'normals_mov'); mesh = rmfield(mesh,'normals_mov'); end;
    if isfield(mesh,'vecprod_mov'); mesh = rmfield(mesh,'vecprod_mov'); end;
    for it = 1:nt
        fprintf('\b\b\b\b\b%02i/%02i',it,nt);
        mesh.qual_mov(:,it) = meshquality(mesh.v_mov(:,:,it),mesh.f);
%        mesh.areas_mov_old(:,it) = elemvolume(mesh.v_mov(:,:,it),mesh.f);
%        mesh.normals_mov_old(:,:,it) = surfacenorm(mesh.v_mov(:,:,it),mesh.f);
%
% --> mgv
        [aaa,nnn,vprod] = compute_vecprod(mesh.v_mov(:,:,it),mesh.f);
      
        mesh.areas_mov(:,it) = aaa;
        if mesh.norm_out
            mesh.normals_mov(:,:,it) = nnn;
            mesh.vecprod_mov(:,:,it) = vprod;
        else
            mesh.normals_mov(:,:,it) = -nnn;
            mesh.vecprod_mov(:,:,it) = -vprod;
        end
    end
    mesh.areas = mesh.areas_mov(:,1);
    fprintf('\nDone.\n'); 
    
elseif nargin==2
    % Use ONE set of moving vertices (v_mov)
    mesh.qual = meshquality(mesh.v_mov(:,:,it),mesh.f);
%    mesh.areas_old = elemvolume(mesh.v_mov(:,:,it),mesh.f);
%    mesh.normals_old = surfacenorm(mesh.v_mov(:,:,it),mesh.f);
    
% --> mgv
        [aaa,nnn,vprod] = compute_vecprod(mesh.v_mov(:,:,it),mesh.f);
        mesh.areas = aaa;
        if mesh.norm_out
            mesh.normals = nnn;
            mesh.vecprod = vprod;
        else
            mesh.normals = -nnn;
            mesh.vecprod = -vprod;
        end

elseif ~isfield(mesh,'v_mov')
    % Use fixed vertices (v)
    mesh.qual = meshquality(mesh.v,mesh.f);
%    mesh.areas_old = elemvolume(mesh.v,mesh.f);
%    mesh.normals_old = surfacenorm(mesh.v,mesh.f);
% --> mgv
        [aaa,nnn,vprod] = compute_vecprod(mesh.v,mesh.f);
        mesh.areas = aaa;
        if mesh.norm_out
            mesh.normals = nnn;
            mesh.vecprod = vprod;
        else
            mesh.normals = -nnn;
            mesh.vecprod = -vprod;
        end

else
    warning(['Trying to compute Mesh Properties for a set of moving '...
             'vertices that does not exist']);
    return;
end
mesh.mean_area = mean(mesh.areas);
return
