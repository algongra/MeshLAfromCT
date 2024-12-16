classdef ODF <  NumberedFile & dynamicprops
%
% \author M.Moriche
% \date 01-05-2015 by M.Moriche \n
%       Created from Frame
%
% \brief Class for HDF5 Output Data Files
%
% \details
%
% Class to simplify the postprocessing of simulation frames.
% Scripts that use this class result much more clear
%
% MANDATORY ARGUMENTS:
%  - basename
%  - path 
%  - ext
%
%
% METHODS:
%  - readfields 
%  - setfield
%  - bytes = getdiskspace
%  
% \code
% path = '/data2/mmoriche/IBcode/...'
% fr = ODF('validation_stcyl');
% fr = ODF('validation_stcyl', path);
% \endcode
%
properties
diskstatus
end

methods

function self = ODF(basenm,path,ext,ifr, varargin)
   % @author M.Moriche
   % @date 01-12-2015 by M.Moriche \n
   %       Created
   %
   % @brief Function to handle output frames from TUCAN
   %
   % @details
   %
   % MANDATORY ARGUMENTS:
   % - basenm: basename of the simulation
   % - path: path where the simulation has been run
   % - ext: extension of the files
   % - ifr': index frame of the simulation to be read. 
   %         Possible values are: integers, 'first', 'last'.
   %
   % OPTIONAL ARGUMENTS:
   % - sizebyvar = struct; default is empty struct.
   % 
   % SETUP:
   % - ifrpos=2: Position of the index in the groups from the SuperClass 
   %             object NumberedFile.regex (Used to read available files)
   %
   % examples:
   %
   %  @code
   %  od = ODF('validation_stcyl','/my/dir','h5frame', 'first')
   %  od = ODF('validation_stcyl','/my/dir','h5frame', 'last')
   %  od = ODF('validation_stcyl','/my/dir','h5frame', 1)
   %  @endcode
   %

   % SuperClass inheritance
   sizebyvar = struct;
   misc.assigndefaults(varargin{:});
   
   ifrpos = 2; % setup for File Object
   self@NumberedFile(path, basenm, ext, ifrpos);

   % get ifr
   r = ODSIM(self.basenm,self.path,self.ext);
   [ifr stat] = r.getifr(ifr);

   self.ifr  = ifr;
   self.diskstatus = stat;
   if stat == 1
      misc.warningmessage('Output Data Field (ODF) not found in disk');
   end

end
   
end  % methods
end  % class

