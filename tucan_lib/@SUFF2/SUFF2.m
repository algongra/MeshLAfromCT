classdef SUFF2 < SUFF 
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
%  - time = steptime
%  - tuni = unitarytime
%  - vorticity 
%  - bytes = getdiskspace
%  
methods
function self = SUFF2(varargin)

   self@SUFF(2, varargin{:});

end

   
end  % methods
end  % class

