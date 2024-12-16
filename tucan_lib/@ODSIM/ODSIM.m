classdef ODSIM < dynamicprops
% \date 01-05-2015 by M.Moriche \n
%       Created

properties
   datapath
   basenm
   nf 
   odfs
   ifrlist
   multidata
end

methods
function self = ODSIM(basenm, datapath, ext, varargin)
   % \date 01-05-2015 by M.Moriche \n
   %       Created
   %function self = ODSIM(basenm, datapath, ext)

   multidata = true;
   misc.assigndefaults(varargin{:});

   self.nf       = NumberedFile(datapath,basenm,ext, 2);
   self.datapath = datapath;
   self.basenm   = basenm;
   self.multidata = multidata;
 
return
end


end %methods

end %classdef
