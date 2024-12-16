classdef FLO223 < SUFF3 
%
% \author M.Moriche
% \date 01-05-2015 by M.Moriche \n
%       Created
properties
   fr2D  % 2D Frame
   frFLO % Floquet Frame
end

methods

function self = FLO223(frFLO, fr2D, varargin);
   
   z0 = 0;
   zf = 1;
   nz = 20;
   misc.assigndefaults(varargin{:});

   self@SUFF3('doslab', false);

   self.fr2D = fr2D;
   self.frFLO = frFLO;

   self.z0 = z0;
   self.zf = zf;
   self.nz = nz;

end

gen3D_z_interp(self, kz);
   
end  % methods
end  % class

