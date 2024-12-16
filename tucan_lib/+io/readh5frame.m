function [fr, failitems] = readh5frame(basenm, varargin)
% @author M.Moriche
% @date Created 06-06-13 by M.Moriche
% @date Modified 07-06-13 by M.Moriche
% @date 07-06-13 by M.Moriche
%       added extra arguments in H5D.read
% @date 22-07-13 by M.Moriche
%       added error handling and failitems as output
% @date 13-12-13 by M.Moriche
%       Documented and organized. Use now readh5file
% 
% @brief Function to read data from a frame saved with HDF5 lib
%
% @details
%
% Uses low level function in HDF5 library. High level libraries are
%  not available in all Matlab releases.
%
% Mandatory arguments: 
%  - basenm: reference name of the simulation
%
% Optional arguments:
%  - ifr: frame index (Default: 'last')
%         can take values 'first', 'last' or integer index of the frame.
%  - items: variables to be read. 
%           Default: {'ux','uy','p','xux','yux','xuy','yuy','xp','yp'}
%  - path where the data is saved.
%
% Examples:
%  @code
%  fr = readh5frame('validation_1')
%  fr = readh5frame('validation_1','ifr',12)
%  fr = readh5frame('validation_1','items',{'ux','uy','phi'})
%  fr = readh5frame('validation_1','path','/data/mmoriche')
%  @endcode
%
% @warning the format of the file is fixed inside this function
%

frmt    = '%s.%05d.h5frame';
regfrmt ='(\w*?)\.(\d+?)\.h5frame';

failitems = {};

ifr = 'last';
items = {'ux','uy','p','xux','yux','xuy','yuy','xp','yp'};
path = '.';
misc.assigndefaults(varargin{:});

ifrclass =  class(ifr); 
if strcmp(ifrclass,'char')
   ifrlist = misc.getifr(fullfile(path,basenm),regfrmt,2);
   if strcmp(ifr,'first')
      ifr = ifrlist(1);
   elseif strcmp(ifr,'last')
      ifr = ifrlist(end);
   end
end

fnm = fullfile(path, basenm,  sprintf(frmt,basenm,ifr));

[fr, failitems] = io.readh5file(fnm,items);

return
end
