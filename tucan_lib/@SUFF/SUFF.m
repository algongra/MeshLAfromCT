classdef SUFF  < handle
% \author M.Moriche
% \date 01-05-2015 by M.Moriche \n
%       Created
%
% \brief Structured Uniform Fluid Field
%
% \details
%
% Class to handle Eulerian Cartesian flow fields.
%
%
% METHODS:
%  - vorticity(dim) sets w-dim vorticity
  
%properties(Access = public)
properties
ndim;
itemsByGroup;
sizebyvar;
slab; 
doslab;

% Primary variables
ux;
uy;
uz;
p;

% Vorticity
wx;
wy;
wz;
Q;

% Eulerian mesh
x0; xf; y0; yf; z0; zf;
nx;ny;nz
dx;dy;dz;
xux;
yux;
zux;
xuy;
yuy;
zuy;
xuz;
yuz;
zuz;
xp;
yp;
zp;

% Boundary conditions
bcp11; bcp12; bcp21; bcp22; bcp31; bcp32;
bcux11; bcux12; bcux21; bcux22; bcux31; bcux32;
bcuy11; bcuy12; bcuy21; bcuy22; bcuy31; bcuy32;
bcuz11; bcuz12; bcuz21; bcuz22; bcuz31; bcuz32;


end

methods

function self = SUFF(ndim,varargin)
% @author M.Moriche
% @date 01-05-2015 by M.Moriche \n
%       Created
%
% @brief Structured Uniform field Class
%
% @details
%
% MANDATORY ARGUMENTS:
% - ndim: rank of the problem (2 or 3)
%
% OPTIONAL ARGUMENTS:
% - slab = NaN(ndim,2); Box to make an slab of the fields
% - doslab = true; flag to Slab the field
% 
% Members:
% 
% 
% examples:
%
% @code
% fr = SUFF(3)
% fr = SUFF(2,'slab',[-1,3;-4,4])
% @endcode
 

% optional arguments
slab = NaN(ndim, 2);
doslab = false;
misc.assigndefaults(varargin{:});

self.ndim      = ndim;
self.slab      = slab;
self.doslab    = doslab;
self.sizebyvar = struct;

% generate items to read fields (dependent on ndim)
itemsByGroup      = self.getgroups(ndim);
self.itemsByGroup = itemsByGroup;
%sizebybar         = self.getsizebyvar();
%self.sizebyvar = sizebyvar;

end

function disp(self)
    tab = {'Structured Uniform Flow Field', ['Dim: ' num2str(self.ndim)]};
    if self.doslab
       row = {'Slabbed', mat2str(self.slab)};
       tab = [tab; row];
    end
    row = {'Groups', misc.strjoin(fields(self.itemsByGroup), ', ')};
    tab = [tab; row];
    mytab.tab2ascii(tab)
return
end

% functions defined in separeted files:
readfields(self,fnm,items)
vorticity(self, dim, varargin) 

end  % methods


end  % class

