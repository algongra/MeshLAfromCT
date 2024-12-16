function [xyz varargout] = old_gmsh2matlab(pthin,flsnmin,h,varargin)
% @author A.Gonzalo
%
% @brief Function to transform ASCII gmsh.msh file/s into matlab variables
%        for create Lagrangian mesh
%
% @date 02-07-2015 by A.Gonzalo \n
%                  Created and documented
%
% @date 02-02-2016 by A.Gonzalo \n
%                  Modified, now it can save connections
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
%  - ifinfo: print and plot information about the mesh. [false by default]
%  - ifplot: plot mesh points colored by associated marker volume.
%            [false by default]
%
% MINIMUM OUTPUT
% --------------
%  - xyz: matrix containing the coordinates of the Lagrangian points
%         whose shape is (n,3), where n is the number of lagrangian points.
%
% OPTIONAL OUTPUT
% ---------------
%  - vol: associated marker volume.
%  - ib(nbody): begin index of each body [integer]
%  - ie(nbody): ending index of each body [integer]
%  - conn: triangles connections (nodes)
%
% EXAMPLES
% --------
%  @code
%  [xyz] = gmsh.gmsh2matlab(pthin,flsnmin,h);
%  [xyz vol] = gmsh.gmsh2matlab(pthin,flsnmin,h);
%  [xyz vol ib ie nbody] = gmsh.gmsh2matlab(pthin,flsnmin,h);
%  [xyz vol ib ie nbody] = gmsh.gmsh2matlab(pthin,flsnmin,h,'ifinfo',true);
%  [xyz vol] = gmsh.gmsh2matlab(pthin,flsnmin,h,'ifplot',true,'ptssize',25);
%  [xyz vol ib ie nbody conn] = gmsh.gmsh2matlab(pthin,flsnmin,h,'ifinfo',true);
%  @endcode

% defaults
ifinfo = false;
ifplot = false;
ptssize = 20;
misc.assigndefaults(varargin{:});

% number of bodies
nbody = length(flsnmin);

i_glb_pts = 0;
n_glb_pts = 0;
etri = 0;
btri = 0;
lmin = 1e22;
lmax = 0;
for ibody = 1:nbody;
    iposconn_rep = [];
    % READ data from gmsh.msh files
    flnm = strcat(flsnmin{ibody},'.msh'); flnm = fullfile(pthin,flnm);
    msh = gmsh.load_gmsh(flnm);
    msh.TRIANGLES = msh.TRIANGLES(:,[1 2 3]);

    % Calculate xyz and vol of all nodes used to build lagrangian mesh and
    % conn (connection between them)
    for inode = 1:msh.nbNod;
        [ipos jpos] = find(msh.TRIANGLES == inode);
        iposconn = setdiff(ipos,iposconn_rep);
        if ~isempty(iposconn)
           % mesh connection
           btri = etri + 1;
           etri = etri + length(iposconn);
           conn(btri:etri,:) = msh.TRIANGLES(iposconn,:);
        end
        if ~isempty(ipos)
           i_glb_pts = i_glb_pts + 1;
           % lagrangian mesh
           xyz(i_glb_pts,1) = msh.POS(inode,1);
           xyz(i_glb_pts,2) = msh.POS(inode,2);
           xyz(i_glb_pts,3) = msh.POS(inode,3);
           conn(conn == inode) = i_glb_pts;
           % marker volume
           [vol(i_glb_pts) mxsz mnsz] = gmsh.get_vol(msh.TRIANGLES(ipos,:),msh.POS,h);
        end
        iposconn_rep = [iposconn_rep ; iposconn];

        % get minimum distance between 2 mesh points
        if lmin > mnsz
           lmin = mnsz;
        end
        % get maximum distance between 2 mesh points
        if lmax < mxsz
           lmax = mxsz;
        end
    end
 
    % begin index of the body
    ib(ibody) = n_glb_pts + 1;
    % ending index of the body
    ie(ibody) = i_glb_pts;
    % number of global points of the mesh
    n_glb_pts = i_glb_pts;

    % clear msh struct and other variables
    clear msh flnm ipos jpos
end

% parameters needed for checks 
maxlpvol = max(vol(:));
eulvol = h^3;
% checks
if maxlpvol > eulvol
   disp(' ')
   disp(' ')
   fprintf('ERROR WRITTING GEOMETRY IN .dat FILE IN %s FUNCTION',mfilename);
   disp(' ')
   disp(' ')
   aggfrmt = ['volume associated to lagrangian points must be smaller than '...
              'eulerian mesh volume'];
   disp(aggfrmt)
   disp(' ')
   fprintf('\nmaximum volume associated to Lagrangian points = %17.15f',maxlpvol);
   fprintf('\neulerian volume (dxr*dyr*dzr)                  = %17.15f',eulvol);
   disp(' ')
   disp(' ')
   quit
end

if ifinfo
   % get minimum value of lagrangian point volume marker
   minlpvol = min(vol(:));
   % print information
   fprintf('\nminimum distance between two Lagrangian points = %17.15f',lmin);
   fprintf('\nmaximum distance between two Lagrangian points = %17.15f',lmax);
   fprintf('\nminimum volume associated to Lagrangian points = %17.15f',minlpvol);
   fprintf('\nmaximum volume associated to Lagrangian points = %17.15f',maxlpvol);
   fprintf('\neulerian volume (dxr*dyr*dzr)                  = %17.15f',eulvol);
   fprintf('\nnumber of lagrangian points                    = %i   \n',n_glb_pts);
   % get lagrangian points volume marker size distribution (in %)
   percvols=vol/maxlpvol*100;
   % plot lagrangian points volume marker size distribution (in %)
   figure(); hist(percvols);
   xlabel('volume/volumen_{max} in %'); ylabel('lagrangian points');
end

if ifplot
   figure(); scatter3(xyz(:,1),xyz(:,2),xyz(:,3),ptssize,vol,'filled');
   axis equal; colorbar;
end

if nargout == 2
       varargout{1} = vol;
elseif nargout == 3
       varargout{1} = vol;
       varargout{2} = conn;
elseif nargout == 5
       varargout{1} = vol;
       varargout{2} = ib;
       varargout{3} = ie;
       varargout{4} = nbody;
elseif nargout == 6
       varargout{1} = vol;
       varargout{2} = ib;
       varargout{3} = ie;
       varargout{4} = nbody;
       varargout{5} = conn;
end

return 
end
