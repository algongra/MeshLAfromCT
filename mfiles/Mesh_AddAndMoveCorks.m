function mesh = Mesh_AddAndMoveCorks(mesh,ncorkplanes,dh,nt_fr)


ncorkite = 40;
n_edgesrem = mesh.n_edges;
edge_nifr = zeros(n_edgesrem-1,3,nt_fr);
fixed_corks = zeros(n_edgesrem-1,nt_fr);

fig = figure();

for it = 1:nt_fr
    for i = 1:n_edgesrem-1
        face_type_indices = (mesh.face_type == i);
        normal_to_inlet = median(mesh.normals_mov(face_type_indices,:,it),1);
        edge_nifr(i,:,it) = normal_to_inlet/norm(normal_to_inlet);
    end
end 
for ied = 1:n_edgesrem-1
    PVfac{ied} = find(mesh.face_type==ied);
end
MVfac = find(mesh.face_type==ied+1);
LAfac = find(mesh.face_type==0);

% Obtain indices of the ring of triangles near the edges of mesh structure faces
% with face_type == 0 (LA wall)
%  - They will be necessary to check if the veins corks are too close to
%    the LA wall points
%  - They will also be necessary to make sure that the displacements of
%    points at the edges are equal to zero
for ied = 1:n_edgesrem
    ifacedge{ied} = [];
    edgerem{ied} = mesh.v_edge_ordered{ied}(:,1);
    for ifac = 1:size(mesh.f(LAfac),1);
        if any(ismember(edgerem{ied},mesh.f(LAfac(ifac),:)))
           ifacedge{ied} = [ifacedge{ied} LAfac(ifac)];
        end
    end
    for it =1:nt_fr
        xyzcifr(ied,:,it) = mean(mesh.v_mov(mesh.v_edge_ordered{ied}(:,1),:,it));  
    end
end

 
% Create PVs corks (checking none of their points are inside the LA wall)
%
for ifr = 1:nt_fr
    disp(sprintf('\n##### Creating corks at Frame %02i #####\n',ifr));

    % Get spatial coordinates of triangles' centers for moved remeshed
    % (repaired) grid with lids (without corks yet)
    clear tric;
    for ifac = 1:length(mesh.f);
        tric(ifac,:) = mean(mesh.v_mov(mesh.f(ifac,:),:,ifr));
    end
    iendc = length(tric);

    % Find indices of points (tric) of the LA wall at a distance
    % smaller than 3*dh from PVs edges points
    for ied = 1:n_edgesrem-1
        iwallclose2edge{ied} = ifacedge{ied}';
        for ic = 1:length(ifacedge{ied})
            dd = sqrt( sum((tric(ifacedge{ied}(ic),:) -...
                            tric(LAfac,:)).^2,2) );
            isin = find( dd < (3*dh) );
            if ~isempty(isin)
               iwallclose2edge{ied} = [iwallclose2edge{ied} ;...
                                       LAfac(isin)];
            end
        end
        iwallclose2edge{ied} = unique(iwallclose2edge{ied});
    end

    % Copy mesh fields (without corks at frame ifr in a temporal mesh
    % (tmpmsh)
    tmpmsh.v = mesh.v_mov(:,:,ifr);
    tmpmsh.qual = mesh.qual_mov(:,ifr);
    tmpmsh.areas = mesh.areas_mov(:,ifr);
    tmpmsh.normals = mesh.normals_mov(:,:,ifr);
    if ifr == 1;
       tmpmsh.f = mesh.f;
       tmpmsh.face_type = mesh.face_type;
    end
    tmpmsh.faux = mesh.f;

    % Add triangulations of corks to tmpmsh (only in the PVs, not in the MV)
    for ipv = 1:n_edgesrem-1
        ifcorktoberepaired = true; ite = 0;

        while ifcorktoberepaired & ite <= ncorkite
          ite = ite + 1;
          if ite == 1
             disp(sprintf('Vein #%i\n Iteration #%i',ipv,ite));
          else
             disp(sprintf(' Iteration #%i',ite));
          end

          % Get last value of faces and vertices indices before including this
          % cork to delete the cork generated in case it was too close to the
          % wall
          iendf = length(tmpmsh.faux); iendv = length(tmpmsh.v);

          % Identify faces for current vein lid
          lidf = mesh.f(PVfac{ipv},:); newf = lidf(:);

          % These quantities will be just copied in each plane of the vein cork
          newa = tmpmsh.areas(PVfac{ipv});
          newq = tmpmsh.qual(PVfac{ipv});
          %newn = tmpmsh.normals(ilidsremf{ipv},:);
          newn = repmat(edge_nifr(ipv,:,ifr),[length(PVfac{ipv}) 1]);

          % Get unique sorted vein lid faces indices
          lidfuni = unique(lidf);

          % Initialize auxiliary structure cork
          cork.v = []; cork.f = [];
          % Add plane to the cork
          for ipl = 1:ncorkplanes
              % Map old vertices onto newly added vertices
              newf = lidf(:);
              for iv = 1:length(lidfuni)
                  newf(newf==lidfuni(iv)) = size(tmpmsh.v,1) + iv;
              end
              newf = reshape(newf,[length(newf)/3 , 3]);

              % Translate vertices a distance dh in the normal direction
              % obtained for the lids (edge_nifr)
              dhvec = repmat(dh,[1 3]);
              newv = bsxfun(@plus,tmpmsh.v(lidfuni,:),...
                                  ipl.*edge_nifr(ipv,:,ifr).*dhvec);

              % Copy newa, newq, newn, and newv in the structure tmpmsh
              tmpmsh.areas   = [tmpmsh.areas ; newa];
              tmpmsh.qual    = [tmpmsh.qual ; newq];
              tmpmsh.normals = [tmpmsh.normals ; newn];
              tmpmsh.v       = [tmpmsh.v ; newv];
              % Copy wf in auxiliar field faux of structure tmpmsh
              tmpmsh.faux    = [tmpmsh.faux ; newf];
              % Copy newf and newv in auxiliary structure cork
              cork.f       = [cork.f ; length(cork.v)+(newf-min(newf)+1)];
              cork.v       = [cork.v ; newv];
              if ifr == 1
                 % Copy newf and newface_type in structure tmpmsh
                 tmpmsh.f         = [tmpmsh.f ; newf];
                 tmpmsh.face_type = [tmpmsh.face_type ;... 
                                     (20+ipv)*ones(size(newf,1),1)];
              end
          end

          disp('  Checking cork...');

          % Find faces indices of the corks
          if ifr == 1
             PVcorkfac{ipv} = find(tmpmsh.face_type==20+ipv);
             corklen(ipv) = length(PVcorkfac{ipv});
          end

          % Get spatial coordinates of triangles' centers for moved remeshed
          % (repaired) grid with lids (with the cork created)
          for ifac = 1:corklen(ipv);
              tric(iendc+ifac,:) = mean(cork.v(cork.f(ifac,:),:));
          end

          % Find indices of points (tric) of the LA wall at a distance smaller
          % than 3*dh from PVs corks points
          iwalltooclose{ipv} = [];
          for ic = 1:corklen(ipv)
              dd = sqrt( sum((tric(iendc+ic,:) -...
                              tric(LAfac,:)).^2,2) );
              isin = find( dd < (3*dh) );
              % Is this point of the cork touching any point of the LA wall not
              % near its associated PV edge?
              if ~isempty(isin) &...
                 ~isempty(setdiff(LAfac(isin),iwallclose2edge{ipv}))
                 isin = setdiff(LAfac(isin),iwallclose2edge{ipv});
                 iwalltooclose{ipv} = [iwalltooclose{ipv} ; isin];
              end
          end
          iwalltooclose{ipv} = unique(iwalltooclose{ipv});
          %
          % Calculate the new value of edge_nifr for this edge at this frame
          if ~isempty(iwalltooclose{ipv})
             %{
             % Plot when a cork is too close to the LA wall
             cc = lines(n_edgesrem);
             set(0,'CurrentFigure',fig); clf; hold on;
             %
             pLAifr = patch('Vertices',tmpmsh.v,'Faces',tmpmsh.f(LAfac,:));
             pLAifr.FaceColor = 0.5*[1 1 1]; pLAifr.FaceAlpha = 0.4;
             view(3); axis equal; daspect([1 1 1]); view(-170,20)
             camlight; lighting gouraud; box on; grid on;
             
             for ii = 1:n_edgesrem-1
                 pPVifr{ii} = patch('Vertices',tmpmsh.v,...
                                    'Faces',tmpmsh.f(PVfac{ii},:));
                 pPVifr{ii}.FaceColor = cc(ii,:);
                 pPVifr{ii}.FaceAlpha = 0.4;
             end
             pMVifr = patch('Vertices',tmpmsh.v,'Faces',tmpmsh.f(MVfac,:));
             pMVifr.FaceColor = cc(n_edgesrem,:);
             %
             p3crk = plot3(tric(iendc+1:end,1),...
                           tric(iendc+1:end,2),...
                           tric(iendc+1:end,3),...
                           '.','Color',cc(ipv,:));
             %
             inwall = patch('Vertices',tmpmsh.v,...
                            'Faces',tmpmsh.f(iwalltooclose{ipv},:),...
                            'FaceColor','red','EdgeColor','none');
             % Plot old normal vector
             ilastcorklayer = round(length(cork.v) - (length(cork.v)/ncorkplanes)/2);
             q3old = quiver3(cork.v(ilastcorklayer,1),...
                             cork.v(ilastcorklayer,2),...
                             cork.v(ilastcorklayer,3),...
                             edge_nifr(ipv,1,ifr),...
                             edge_nifr(ipv,2,ifr),...
                             edge_nifr(ipv,3,ifr),...
                             2,'Color','r','LineWidth',2);
             pause(0.01);
             %}
             disp('   Too close to the LA wall');
             %
             % Fit a plane to the points of the LA wall that are too close to
             % the cork (iwalltooclose{ipv})
             if length(iwalltooclose{ipv}) > 2
                [npl,~,ppl] = affine_fit(tric(iwalltooclose{ipv},:));
                npl = npl';
                if dot(npl,mean(tmpmsh.normals(iwalltooclose{ipv},:)))<0;
                   npl = -npl;
                end
             elseif length(iwalltooclose{ipv}) == 2
                npl = mean(tmpmsh.normals(iwalltooclose{ipv},:));
                ppl = mean(tric(iwalltooclose{ipv},:));%??????????
             else
                npl = tmpmsh.normals(iwalltooclose{ipv},:);
                ppl = tric(iwalltooclose{ipv},:);
             end
             %{
             % Plot normal vector perpendicular to the points
             q3per = quiver3(ppl(1),ppl(2),ppl(3),npl(1),npl(2),npl(3),...
                             2,'Color','b','LineWidth',2);
             %}
             %
             nplrep = repmat(npl,[length(tric(iwalltooclose{ipv})) 1]);
             %{
             distn = max(abs(dot(tric(iwalltooclose{ipv},:)-ppl,nplrep,2))) +...
                     dh;
             %}
             distxyzc2ppl = sqrt(sum((ppl-xyzcifr(ipv,:,ifr)).^2));
             dd = dot(tric(iwalltooclose{ipv},:) - xyzcifr(ipv,:,it),nplrep,2);
             distn = max(dd) - min(dd);
             % When only one point is detected in iwalltooclose{ipv}, distn is 
             % equal to zero and naux is not modified (updated)
             % Similarly, when the cloud of points in iwalltooclose{ipv} is
             % small, the algorithm cannot rotate the cork in ncorkite
             if exist('distnold')
                if distn == distnold
                   distn = distn + dh/4;
                end
             end
             distnold = distn;
             %
             naux = distxyzc2ppl*edge_nifr(ipv,:,ifr) + distn*npl;
             naux = naux/norm(naux);
             % Overwrite edge_nifr
             edge_nifr(ipv,:,ifr) = naux;
             fixed_corks(ipv,ifr) = 1;
             %{
             % Plot new normal vector
             q3new = quiver3(cork.v(ilastcorklayer,1),...
                             cork.v(ilastcorklayer,2),...
                             cork.v(ilastcorklayer,3),...
                             edge_nifr(ipv,1,ifr),...
                             edge_nifr(ipv,2,ifr),...
                             edge_nifr(ipv,3,ifr),...
                             2,'Color','g','LineWidth',2);
             pause(0.01);
             %}
             % Delete cork faces and edges created
             tmpmsh.areas(iendf+1:end)     = [];
             tmpmsh.qual(iendf+1:end)      = [];
             tmpmsh.normals(iendf+1:end,:) = [];
             tmpmsh.v(iendv+1:end,:)       = [];
             if ifr == 1
                tmpmsh.f(iendf+1:end,:)       = [];
                tmpmsh.face_type(iendf+1:end) = [];
             end
             tmpmsh.faux(iendf+1:end,:)       = [];
          else
             %{
             % Plot correct cork
             cc = lines(n_edgesrem);
             set(0,'CurrentFigure',fig); clf; hold on;
             %
             pLAifr = patch('Vertices',tmpmsh.v,'Faces',tmpmsh.f(LAfac,:));
             pLAifr.FaceColor = 0.5*[1 1 1]; pLAifr.FaceAlpha = 0.4;
             view(3); axis equal; daspect([1 1 1]);
             camlight; lighting gouraud; box on; grid on;
             
             for ii = 1:n_edgesrem-1
                 pPVifr{ii} = patch('Vertices',tmpmsh.v,...
                                    'Faces',tmpmsh.f(PVfac{ii},:));
                 pPVifr{ii}.FaceColor = cc(ii,:);
                 pPVifr{ii}.FaceAlpha = 0.4;
             end
             %
             pMVifr = patch('Vertices',tmpmsh.v,'Faces',tmpmsh.f(MVfac,:));
             pMVifr.FaceColor = cc(n_edgesrem,:);
             pMVifr.FaceAlpha = 0.4;
             %
             corkfac{ipv} = find(tmpmsh.face_type==20+ipv);
             pcrk = patch('Vertices',tmpmsh.v,'Faces',tmpmsh.f(corkfac{ipv},:));
             pcrk.FaceColor = cc(ipv,:);
             pcrk.FaceAlpha = 0.4;
             pause(0.01);
             %}
             disp('   Valid');
             %
             ifcorktoberepaired = false;
          end
          % Delete tric faces of the cork created
          tric(iendc+1:end,:) = [];
        end
        if ite > ncorkite
           errortxt = ['\nERROR: Displacements values from original grids '...
                       'have not been read correctly\n'];
           error(sprintf(errortxt));
        end
    end
    % Copy moving mesh with corks of frame ifr on *mov fields of tmpmsh
    tmpmsh.areas_mov(:,ifr) = tmpmsh.areas;
    tmpmsh.qual_mov(:,ifr) = tmpmsh.qual;
    tmpmsh.normals_mov(:,:,ifr) = tmpmsh.normals ;
    tmpmsh.v_mov(:,:,ifr) = tmpmsh.v;


    if ifr == 1
       for ied = 1:n_edgesrem-1
           PVcorkfac{ied} = find(tmpmsh.face_type==20+ied);
       end
    end

    %{
    % Plot LA wall at this frame
    set(0,'CurrentFigure',fig); clf; hold on;
    pLArem = patch('Vertices',tmpmsh.v_mov(:,:,ifr),...
                   'Faces',tmpmsh.f(LAfac,:));
    pLArem.FaceColor = 0.25*[1 1 1]; pLArem.FaceAlpha = 1;
    view(3); axis equal; daspect([1 1 1]); camlight; lighting gouraud;
    box on; grid on;
    % Plot PVs, PVs corks, and MV at this frame
    for ipv = 1:n_edgesrem-1
        ppvlid{ipv} = patch('Vertices',tmpmsh.v_mov(:,:,ifr),...
                            'Faces',tmpmsh.f(PVfac{ipv},:));
        ppvlid{ipv}.FaceColor = cc(ipv,:);
        ppvlid{ipv}.FaceAlpha = 0.4;
        %
        ppvcrk{ipv} = patch('Vertices',tmpmsh.v_mov(:,:,ifr),...
                            'Faces',tmpmsh.f(PVcorkfac{ipv},:));
        ppvcrk{ipv}.FaceColor = cc(ipv,:);
        ppvcrk{ipv}.FaceAlpha = 0.4;
    end
    pmvlid = patch('Vertices',tmpmsh.v_mov(:,:,ifr),'Faces',tmpmsh.f(MVfac,:));
    pmvlid.FaceColor = cc(n_edgesrem,:); 
    pmvlid.FaceAlpha = 0.4;
    % Plot edges of PVs and MV
    for ied = 1:n_edgesrem
        p3ed{ied} = plot3(tmpmsh.v_mov(edgerem{ied},1,ifr),...
                          tmpmsh.v_mov(edgerem{ied},2,ifr),...
                          tmpmsh.v_mov(edgerem{ied},3,ifr),...
                          'LineWidth',2,'Color',cc(ied,:));
    end
     Plot PVs normal vectors at this frame
    for ipv = 1:n_edgesrem-1
        xyzcaux = xyzcifr(ipv,:,ifr);
        edgenaux = edge_nifr(ipv,:,ifr);
        qnedpl{ipv} = quiver3(xyzcaux(1),xyzcaux(2),xyzcaux(3),...
                              edgenaux(1),edgenaux(2),edgenaux(3),...
                              'LineWidth',2,'Color',cc(ipv,:));
    end
    %}
end
% Delete *mov fields from mesh structure to copy the new values
mesh = rmfield(mesh,'areas_mov');
mesh = rmfield(mesh,'normals_mov');
mesh = rmfield(mesh,'qual_mov');
mesh = rmfield(mesh,'v_mov');
mesh = rmfield(mesh,'f');
mesh = rmfield(mesh,'face_type');
% Copy new values of *mov fields
mesh.areas_mov   = tmpmsh.areas_mov;
mesh.qual_mov    = tmpmsh.qual_mov;
mesh.normals_mov = tmpmsh.normals_mov;
mesh.v_mov       = tmpmsh.v_mov;
mesh.f           = tmpmsh.f;
mesh.face_type   = tmpmsh.face_type;
% Create new variable to label edges' in which corks have been fixed
mesh.fixed_corks = fixed_corks;
% Create new variable for edges' normals (avoiding corks to pierce LA surface)
mesh.normals_edges = edge_nifr;

return

end

