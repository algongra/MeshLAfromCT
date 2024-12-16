function [ifrlist stat]= getifrlist(path,regfrmt,pos)
% @author M.Moriche
% @date 13-12-2013 by M.Moriche \n
%       Created
%
% @brief Function to read files in path $path saved with 
%        regular expression format $regfrmt where the index
%        of the frame is the group pos of the regular expression.
%
% @details
%
% Lists content of path and checks every file inside to 
% match the regular expression $regfrmt. If it matches, 
% extracts the token at position $pos and converts it to a
% number. This number must define the index of the frame.
%
% Mandatory arguments
%   - path: path to list the content.
%   - regfrmt: regular expression with the appropiate tokens.
%   - post: token index of the frame in the regular expression.
%
% NO Optional arguments
%
% @warning To fast the routine and avoid reallocations, a 
%          maximum size of ifrmax is defined inside to allcoate
%          ifrlist. Then ifrlist is cut to the proper size.
%

listing = dir(path);

n = length(listing);


if n == 0
   % if the directory is not found or it is empty
   stat = 1;
   ifrlist = [];
else
   % if the directory exists
   ifrmax = 10000;
   ifrlist = zeros(1,ifrmax);
   
   i2=1;
   for i = 1:n
      fnm = listing(i).name;
      a = regexp(fnm,regfrmt); 
      if a
         a = regexp(fnm,regfrmt,'tokens'); 
         ifr = a{1}{pos};
         ifrlist(i2) = str2num(ifr);
         i2 = i2 + 1;
      end
   end
   
   if i2 == 1;
      % if  there were no mathcings
      stat = 1;
      ifrlist = [];
   else
      stat = 0;
      ifrlist = ifrlist(1:i2-1);
   end
end

return
end
