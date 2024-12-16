function [commonvalues commonpos] = getcommon(varargin)
% @author M.Moriche
% 
% @brief Function to get common items in serval arrays
%
% @date 26-12-2013 by M.Moriche \n
%       Created
% 
% @details 
%
% Finds the items that sevar arrays have in common. 
%
% OUTPUT
% - commonvalues: the common items
% - commonpos: cell array giving the positions of the common
%              items for each
%
% MANDATORY ARGUMENTS
% - at least two arrays should be input
%
% NO OPTIONAL ARGUMENTS
%
%
% EXAMPLES:
%
% [comm, pos] = getcommon(array1, array2, array3)
% [comm, pos] = getcommon(array1, array2, array3, array4)
%

commonvalues = [];
commonpos = [];

for i = 1:nargin
   arr1 = varargin{i};
   ifound = 1;
   for k = 1:length(arr1)
      arr1(k);
      inall = true;
      counter = 0;
      for j = 1:nargin 
          if i ~= j
             arr2 = varargin{j};
             uu = find(arr1(k) == arr2);
             if ~isempty(uu)
                counter = counter + 1;
             end
          end
      end
      if counter == nargin-1;
         if isempty(find(commonvalues == arr1(k)))
            commonvalues = [commonvalues arr1(k)];
         end
         commonpos(i,ifound) = k;
         ifound = ifound + 1;
      end
   end
end


return
end
