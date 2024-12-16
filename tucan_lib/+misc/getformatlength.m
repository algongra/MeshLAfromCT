function [narray varargout] = getformatlength(patt, varargin)

% OPTIONAL ARGUMENTS

%  path
path = '.';
misc.assigndefaults(varargin{:});

filelist = dir(path);

narray = 0;
stat = 1;

for i = 1:length(filelist)
   f = filelist(i);
   fnm = f.name;
   a = regexp(fnm, patt, 'start');
   if ~isempty(a)
      stat = 0;
      b = regexp(fnm, patt, 'tokens');
      narray = zeros(1,length(b{1})); 
      for j = 1:length(b{1})
        narray(j) = length(b{1}{j});
      end
      break
   end
end

varargout{1} = stat;

return
end
