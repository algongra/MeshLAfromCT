function [vol varargout] = get_vol(tri,gn,h)
% @author A.Gonzalo
%
% @brief Function to calculate the volume of each gmsh node.
%
% @date 08-05-2015 by A.Gonzalo \n
%                  Created
% @date 08-05-2015 by A.Gonzalo \n
%                  Documented
%
% @details
%
% MANDATORY ARGUMENTS
% -------------------
%  - tri: vertices of triangles created by gmsh node. [?x3 double]
%         ? is the number of triangles created by the gmsh node selected
%         This array contains the connectivity (indices of the gmsh nodes
%         of the vertices of the triangles).
%  - gn: coordinates of gmsh nodes. [?x3 double]
%        ? is the number of gmsh nodes of the mesh.
%  - h: minimum space between eulerian mesh points. [double]
%
% MINIMUM OUTPUT
% --------------
%  - vol: associated marker volume.
%
% OPTIONAL OUTPUT
% ---------------
%  - maxsize: maximum distance between 2 gmsh nodes.
%  - minsize: minimum distance between 2 gmsh nodes.
%
% EXAMPLE
% -------
%  @verbatim
%  % Read data from mygmshmesh.msh file located in /home/user/gmshmeshes
%  msh = gmsh.load_gmsh('/home/user/gmshmeshes/mygmshmesh.msh');
%
%  % Get an array containing the indices of the nodes of the vertices of the
%  % triangles
%  msh.TRIANGLES = msh.TRIANGLES(:,[1 2 3]);
%
%  % Get the vertices of the triangles built by the gmsh node 12
%  inode = 12;
%  [ipos jpos] = find(msh.TRIANGLES == 12);
%  vertices = mesh.TRIANGLES(ipos,:);
%
%  % USAGE OF THIS FUNCTION
%  h = 0.01171875;
%  [vol(inode)] = gmsh.get_vol(vertices,mesh.POS,h);
%  [vol(inode) mxsz mnsz] = gmsh.get_vol(vertices,mesh.POS,h);
%  @endverbatim

gn_surf = 0;
ntri = size(tri,1);
minsize = 1e22;
maxsize = 0;
for itri = 1:ntri;
    % get the coordinates of the triangle vertices
    xA = gn(tri(itri,1),1); yA = gn(tri(itri,1),2); zA = gn(tri(itri,1),3);
    xB = gn(tri(itri,2),1); yB = gn(tri(itri,2),2); zB = gn(tri(itri,2),3);
    xC = gn(tri(itri,3),1); yC = gn(tri(itri,3),2); zC = gn(tri(itri,3),3);

    % calculate the distance of the triangle sides
    a = sqrt((xB-xC)^2 + (yB-yC)^2 + (zB-zC)^2);
    b = sqrt((xA-xC)^2 + (yA-yC)^2 + (zA-zC)^2);
    c = sqrt((xA-xB)^2 + (yA-yB)^2 + (zA-zB)^2);

    % calculate triangle surface with Heron's formula
    s = (a+b+c)/2;
    tri_surf = sqrt(s*(s-a)*(s-b)*(s-c));

    % calculate gmsh node associated surface
    gn_surf = gn_surf + tri_surf/3;

    % get minimum distance between 2 gmsh nodes
    dmmy = min(a,b); dmmy = min(dmmy,c);
    if minsize > dmmy
       minsize = dmmy;
    end
    % get maximum distance between 2 gmsh nodes
    dmmy = max(a,b); dmmy = max(dmmy,c);
    if maxsize < dmmy
       maxsize = dmmy;
    end
end
vol = gn_surf*h;

if nargout == 3
   varargout{1} = maxsize;
   varargout{2} = minsize;
end
