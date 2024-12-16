classdef SNUFF3 < SUFF3
%
% \author A.Gonzalo
% \date 06-05-2016 by A.Gonzalo \n
%       Created from SUFF3
%
% \brief Structured Non-Uniform Fluid Field
%
% \details
%
% Class to simplify the postprocessing of simulation frames.
% Scripts that use this class result much more clear
%
% METHODS:
%  - vorticity 
%  - getQ
%  - getgroups
%

%properties(Access = public)
properties
nps;

% Primary variables
ps;

% Eulerian mesh
dxr;
dyr;
dzr;

% Boundary conditions
bcps11; bcps12; bcps21; bcps22; bcps31; bcps32;


end  % properties
  
methods
function self = SNUFF3(varargin)
% @author A.Gonzalo
% \date 06-05-2016 by A.Gonzalo \n
%       Created from SUFF3
%
% @brief Structured Non-Uniform field Class
%
% @details
%
% MANDATORY ARGUMENTS:
%                                     
%
% OPTIONAL ARGUMENTS:
% - slab = NaN(ndim,2); Box to make an slab of the fields
% - doslab = true; flag to Slab the field
% - nps = 0; number of passive scalars fields solved
% 
% Members:
% 
% 
% examples:
%
% @code
% fr = SNUFF()
% fr = SNUFF('nps',1)
% @endcode


% SuperClass inheritance
self@SUFF3(varargin{:});

% Optionals arguments
nps = 0;
misc.assigndefaults(varargin{:});

self.nps = nps;


% generate items to read fields (dependent on ndim)
itemsByGroup      = self.getgroups(self.ndim);
self.itemsByGroup = itemsByGroup;

end

function disp(self)
    tab = {'Structured Non-Uniform Flow Field', ['Dim: ' num2str(self.ndim)]};
    if self.doslab
       row = {'Slabbed', mat2str(self.slab)};
       tab = [tab; row];
    end
    if self.nps > 0
       row = {'Number of passive scalars fields', mat2str(self.nps)};
       tab = [tab; row];
    end
    row = {'Groups', misc.strjoin(fields(self.itemsByGroup), ', ')};
    tab = [tab; row];
    mytab.tab2ascii(tab)
return
end

vorticity(self, idir)
getQ(self, varargin)
   
end  % methods
end  % class

