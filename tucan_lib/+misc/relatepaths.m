function [relativepath isrel] = relatepaths(fullpath, shortpath)
%function [relativepath isrel] = relatepaths(fullpath, shortpath)
% @author M.Moriche 
% @date 29-10-2014
% @brief Function to relate to paths and give the relative one
%
% @details 
%
% OUTPUT: 
%  - relativepath: relaative path from short path to fullpath
%  - isrel: boolean indicating if short path is contained in fullpath
%
% MANDATORY ARGUMENTS: 
%  - fullpath: long path
%  - shortpath: long path
%
%
% @warning If isrel == false, the results is unpredictible
%
%
% @code
% objectspath = '/data2/mmoriche/myproject/notes/myfigures/example1';
% notespath   = '/data2/mmoriche/myproject/notes';
% [relativepath isrel] = relatepaths(objectspath, notespath);
% @endcode
%


aa = misc.strsplit(fullpath,  filesep);
bb = misc.strsplit(shortpath, filesep);

cc = {};
n2 = 0;
for i1 = 1:length(bb)
   if strcmp(aa{i1},bb{i1})
      n2 = n2 + 1;
   end
end

for i1 = length(bb)+1:length(aa)
   cc = [cc aa{i1}];
end

n = length(bb);
isrel = true;
if ~ (n2 == n) 
   isrel = false;
end

relativepath = fullfile(cc{:});

return
end
