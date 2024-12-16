function [myminpos mymaxpos] = findminmax(myarray)
% @author M.Moriche
%
% @brief Function to get the positions of relative maximums and 
%        minimums of an array.
%
% @date 26-12-2013 by M.Moriche \n
%       Taken from old function and documented.
%
% @details
%
% No signal processing, which means that every time the array
% has a change in his growing behavior, the functions detects it
% as a relative maximum/minimum.
%
% OUTPUT
% - myminpos: position in the array of the relative minimums.
% - mymaxpos: position in the array of the relative maximums.
%
% MANDATORY ARGUMENTS
% - myarray: array to look for relative maximum, minimum.
%
% EXAMPLE:
%
% [imin imax] = findminmax(myarray)
%

mymin = [];mymax = [];
%
mydiff=diff(myarray);
myprod=mydiff(1:end-1).*mydiff(2:end);
myvalues1 = myprod<=0;
myvalues2 = mydiff(1:end-1) < 0;
myvalues = myvalues1.*myvalues2;
myminpos=find(myvalues)+1;
%
myvalues2 = mydiff(1:end-1) > 0;
myvalues = myvalues1.*myvalues2;
mymaxpos=find(myvalues)+1;
%
return
end

