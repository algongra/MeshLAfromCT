function [xyz varargout] = gmsh2matlab(pthin,flsnmin,h,varargin)
% @author A.Gonzalo
%
% @brief Function to transform ASCII gmsh.msh file/s into matlab variables
%        for create lagrangian mesh
%
% @date 02-07-2015 by A.Gonzalo \n
%                  Created and documented
%
% @date 02-02-2016 by A.Gonzalo \n
%                  Modified
%                  Connectivity between lagrangian points is saved.
%
% @date 13-07-2016 by A.Gonzalo \n
%                  Modified
%                  lagrangian points are located at the middle of triangles,
%                  instead of at triangles' vertices.
%                  Coordinates of triangles' vertices and its connectivity are
%                  saved (optionally) in xyzv and conn.
%                  Normal vectors associated to the lagrangian points are saved
%                  (optionally) in normvec.
%
% @date 18-08-2016 by A.Gonzalo \n
%                  Modified
%                  lagrangian points can be generated in specific surfaces of
%                  flsnmin through the optional argument surf_type.
%
% @date 30-09-2016 by A.Gonzalo \n
%                  Modified
%                  Added meshflag to check feasibility of the mesh
%
% @details
%
% - Each gmsh.msh file must contain only one body (for create several bodies,
%   several gmsh.msh files must be created).
%
% MANDATORY ARGUMENTS
% -------------------
%  - pthin: path where gmsh.msh files are saved. [string]
%  - flsnmin: name of gmsh.msh files (without .msh). [cell with strings]
%  - h: minimum space between eulerian mesh points. [double]
%
% OPTIONAL ARGUMENTS
% ------------------
%  - ifinfo: print and plot information about the mesh.
%            Size: 1x1.
%            Class: logical.
%            Default: false.
%  - ifplot: plot mesh and colored it as function of the marker volumes
%            associated to each lagrangain point.
%            Size: 1x1.
%            Class: logical.
%            Default: false.
%  - xyzn: point to calculate direction of normal vectors.
%          Size: nbodyx3.
%          Class: double.
%          Default: zeros(nbody,3).
%  - surf_type: flag to select surfaces of each gmsh.msh file where lagrangian
%               points are going to be generated [By default it will get points
%               from all surfaces of flsnmin].
%               Size: 1xnbody.
%               Class: cell.
%               Default: cell(1,nbody).
%
% MINIMUM OUTPUT
% --------------
%  - xyz: coordinates of lagrangian points.
%         Size: nx3. Where n is the number of lagrangian points.
%         Class: double.
%
% OPTIONAL OUTPUT
% ---------------
%  - vol: marker volume associated to each lagrangian point.
%         Size: 1xn. Where n is the number of lagrangian points.
%         Class: double.
%  - ib: begin index of each mesh's body.
%        Size: 1xnbody.
%        Class: integer.
%  - ie: ending index of each mesh's body.
%        Size: 1xnbody.
%        Class: integer.
%  - nbody: number of bodies of the lagrangian mesh.
%           Size: 1x1.
%           Class: integer.
%  - normvec: normal vector associated to each lagrangian point.
%             Size: nx3. Where n is the number of lagrangian points.
%             Class: double.
%  - xyzv: coordinates of triangles's vertices.
%          Size: nvx3. Where nv is the number of triangles' vertices.
%          Class: double.
%  - conn: triangles connections (indices of each triangles's vertices).
%          Size: nx3. Where n is the number of lagrangian points.
%          Class: double.
%  - meshflag: flag to specify if lagrangian mesh fulfill the requirements of
%              a TUCAN's valid mesh.
%              [True. When marker volumes of all lagrangian point are smaller
%              or equal than the Eulerian mesh volume].
%              Size: 1x1.
%              Class: logical.
%
% EXAMPLES
% --------
%  @code
%  flsnmin = {mesh1,mesh2}; pthin = './path_lvl1/path_lvl2'; h = 1/56;
%  [xyz] = gmsh.gmsh2matlab(pthin,flsnmin,h);
%  [xyz vol] = gmsh.gmsh2matlab(pthin,flsnmin,h);
%  xyzn = [0.5 0 2.25];
%  [xyz vol normvec] = gmsh.gmsh2matlab(pthin,flsnmin,h,'xyzn',xyzn);
%  [xyz vol ib ie nbody] = gmsh.gmsh2matlab(pthin,flsnmin,h);
%  [xyz vol ib ie nbody normvec] = gmsh.gmsh2matlab(pthin,flsnmin,h);
%  [xyz vol ib ie nbody normvec xyzv conn] = gmsh.gmsh2matlab(pthin,flsnmin,h);
%  [xyz vol] = gmsh.gmsh2matlab(pthin,flsnmin,h,'ifinfo',true,'ifplot',true);
%  surf_type = {[1 2] [1 3 4]};
%  [xyz vol normvec] = gmsh.gmsh2matlab(pthin,flsnmin,h,'surf_type',surf_type);
%  @endcode

% number of bodies
nbody = length(flsnmin);

% defaults
ifinfo = false;
ifplot = false;
xyzn = zeros(nbody,3);
surf_type = cell(1,nbody);
misc.assigndefaults(varargin{:});

meshflag = false;

% lagrangian mesh required initializations
i_glb_pts = 0;
n_glb_pts = 0;
xver = zeros(1,3); yver = zeros(1,3); zver = zeros(1,3);
lsid = zeros(1,3);
vt = zeros(1,3);
lmin = realmax;
lmax = realmin;


for ibody = 1:nbody;
    % READ data from gmsh.msh files
    %
    flnm = strcat(flsnmin{ibody},'.msh'); flnm = fullfile(pthin,flnm);
    msh = gmsh.load_gmsh(flnm);
    msh.TRIANGLES = msh.TRIANGLES(:,[1 2 3]);
    msh.TRIANGLES_TAGS = msh.ELE_TAGS(msh.ELE_INFOS(:,2)==2,:);

    if ~isempty(surf_type{ibody})
       dmmy = [];
       for itype=1:length(surf_type{ibody})
           dmmy = [dmmy ; msh.TRIANGLES(msh.TRIANGLES_TAGS==surf_type{ibody}(itype),:)];
       end
       msh.TRIANGLES = dmmy;
       clear dmmy
    end

    % Calculate xyz, vol and normvec of lagrangian points (located in the middle
    % of triangles)
    %
    for itri = 1:length(msh.TRIANGLES)
        i_glb_pts = i_glb_pts + 1;

        % coordinates of triangle's vertices
        %
        for isid = 1:3
            xver(isid) = msh.POS(msh.TRIANGLES(itri,isid),1);
            yver(isid) = msh.POS(msh.TRIANGLES(itri,isid),2);
            zver(isid) = msh.POS(msh.TRIANGLES(itri,isid),3);
        end

        % lagrangian point coordinates (xyz)
        %
        xyz(i_glb_pts,1) = mean(xver);
        xyz(i_glb_pts,2) = mean(yver);
        xyz(i_glb_pts,3) = mean(zver);

        % length of triangle's sides
        %
        lsid(1) = sqrt( (xver(2)-xver(3))^2 +...
                        (yver(2)-yver(3))^2 +...
                        (zver(2)-zver(3))^2 );
        lsid(2) = sqrt( (xver(1)-xver(3))^2 +...
                        (yver(1)-yver(3))^2 +...
                        (zver(1)-zver(3))^2 );
        lsid(3) = sqrt( (xver(1)-xver(2))^2 +...
                        (yver(1)-yver(2))^2 +...
                        (zver(1)-zver(2))^2 );

        % marker volume (vol)
        %
        s = sum(lsid)/2;
        %
        vol(i_glb_pts) = sqrt( s*(s-lsid(1))*(s-lsid(2))*(s-lsid(3)) )*h;

        % normal vector (normvec)
        %
        v1 = [xver(2)-xver(1), yver(2)-yver(1), zver(2)-zver(1)];
        v2 = [xver(1)-xver(3), yver(1)-yver(3), zver(1)-zver(3)];
        %
        v1v2 = cross(v1,v2);
        normvec(i_glb_pts,:) = v1v2/norm(v1v2);
        %
        vt = [xyz(i_glb_pts,1)-xyzn(ibody,1),...
              xyz(i_glb_pts,2)-xyzn(ibody,2),...
              xyz(i_glb_pts,3)-xyzn(ibody,3)];
        %
        if dot(vt,normvec(i_glb_pts,:)) < 0
           normvec(i_glb_pts,:) = -normvec(i_glb_pts,:);
        end
    end


    % begin index of the body
    %
    ib(ibody) = n_glb_pts + 1;


    % ending index of the body
    %
    ie(ibody) = i_glb_pts;


    % number of global points of the mesh
    %
    n_glb_pts = i_glb_pts;


    if nargout == 9 | ifplot
       % connectivity of triangles (ndices of each triangles's vertices)
       %
       conn = msh.TRIANGLES;

       % coordinates of triangle's vertices
       %
       xyzv = msh.POS;
    end


    % clear msh struct and name of gmsh.msh read (just in case)
    %
    clear msh flnm
end


% parameters required to assess the feasibility of the mesh
%
maxlpvol = max(vol(:));
eulvol = h^3;


% check feasibility of the mesh
%
if maxlpvol > eulvol

   dmmyfrmt = '\nERROR WRITTING GEOMETRY IN .dat FILE IN %s FUNCTION';
   disp(sprintf(dmmyfrmt,mfilename));
   dmmyfrmt = ['\nVolumes associated to Lagrangian points must be smaller '...
               'than Eulerian mesh volume'];
   disp(sprintf(dmmyfrmt));
   dmmyfrmt = '\nMaximum volume associated to a Lagrangian point = %g';
   disp(sprintf(dmmyfrmt,maxlpvol));
   dmmyfrmt = 'Eulerian mesh volume (dxr*dyr*dzr)              = %g\n';
   disp(sprintf(dmmyfrmt,eulvol));

else
   meshflag = true

   if ifinfo
      % get minimum value of lagrangian point volume marker
      %
      minlpvol = min(vol(:));
   
      % print information
      %
      dmmyfrmt = '\nMinimum volume associated to a Lagrangian point = %g';
      disp(sprintf(dmmyfrmt,minlpvol));
      dmmyfrmt = 'Maximum volume associated to a Lagrangian point = %g';
      disp(sprintf(dmmyfrmt,maxlpvol));
      dmmyfrmt = 'Eulerian mesh volume (dxr*dyr*dzr)              = %g';
      disp(sprintf(dmmyfrmt,eulvol));
      dmmyfrmt = '\nNumber of Lagrangian points                     = %i\n';
      disp(sprintf(dmmyfrmt,n_glb_pts));
   
      % get lagrangian points volume marker size distribution (in %)
      %
      percvols=vol/maxlpvol*100;
   
      % plot lagrangian points volume marker size distribution (in %)
      %
      figinfo1 = figure();
      set(0,'CurrentFigure',figinfo1);
      hist(vol); hold on;
      dmmy = ylim;
      plot(eulvol.*[1 1],dmmy,'-ko');
      hold off;
      xlabel('volume'); ylabel('lagrangian points');

      figinfo2 = figure();
      set(0,'CurrentFigure',figinfo2);
      hist(percvols);
      xlabel('volume/volumen_{max} in %'); ylabel('lagrangian points');
   end
   
   if ifplot
      figmesh = figure();
      set(0,'CurrentFigure',figmesh);
      hold on; view(3);
      patch('Vertices',xyzv,'Faces',conn,'FaceVertexCData',vol','FaceColor','flat')
      axis equal; colorbar;
      xlabel('x/c'); ylabel('y/c'); zlabel('z/c');
   end
   
end

% assign variables arguments go the output
%
if nargout == 2
       varargout{1} = meshflag;
elseif nargout == 3
       varargout{1} = vol;
       varargout{2} = meshflag;
elseif nargout == 4
       varargout{1} = vol;
       varargout{2} = normvec;
       varargout{3} = meshflag;
elseif nargout == 6
       varargout{1} = vol;
       varargout{2} = ib;
       varargout{3} = ie;
       varargout{4} = nbody;
       varargout{5} = meshflag;
elseif nargout == 7
       varargout{1} = vol;
       varargout{2} = ib;
       varargout{3} = ie;
       varargout{4} = nbody;
       varargout{5} = normvec;
       varargout{6} = meshflag;
elseif nargout == 9
       varargout{1} = vol;
       varargout{2} = ib;
       varargout{3} = ie;
       varargout{4} = nbody;
       varargout{5} = normvec;
       varargout{6} = xyzv;
       varargout{7} = conn;
       varargout{8} = meshflag;
end

return 
end
