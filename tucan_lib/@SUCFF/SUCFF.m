classdef SUCFF < SUFF 
%
% \author M.Moriche
% \date 01-05-2015 by M.Moriche \n
%       Created
%
% \brief Structured Uniform Fluid Field
%
% \details
%
% Class to simplify the postprocessing of simulation frames.
% Scripts that use this class result much more clear
%
% METHODS:
%  - C2R  
%

properties
ux_re; ux_im;
uy_re; uy_im;
uz_re; uz_im;
 p_re;  p_im;

end

methods

function self = SUCFF(ndim, varargin)

   % SuperClass inheritance
   self@SUFF(ndim, varargin{:});

   % generate items to read fields (dependent on ndim)
   itemsByGroup2 = self.getcmplxgroups(ndim);

   fldlist = fields(itemsByGroup2);
   for i0 = 1:length(fldlist)
      fld = fldlist{i0};
      self.itemsByGroup.(fld) = itemsByGroup2.(fld);
   end

end

C2R_z( beta, z, varargin);
   
end  % methods
end  % class

