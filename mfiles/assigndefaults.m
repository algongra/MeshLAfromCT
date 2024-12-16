function assigndefaults(varargin)
% @author M.Moriche 
% @date 25-05-2013 by M.Moriche \n
%       Created
% @date 17-12-2013 by M.Moriche \n
%       Documented
%
% @brief Assign default arguments to functions
%
% @details
%
% Assigns value to the caller function workspace.
%
% The syntax of the arguments must be
% @code
% varargin = {variable1name,value1,variable2name,value2}
% @endcode
%
% The names, variable1 and variable2 must be character strings
%
% Example:
%
% @code
% function myfunction(arg1,arg2, varargin)
% 
% x2 = 20;
% alpha = 1;
%
% assigndefaults(varargin{:})
%
% end
% @endcode
% 
% Called by
%
% @code
% myfunction(arg1,arg2)
% myfunction(arg1,arg2,'x2',10)
% myfunction(arg1,arg2,'x2',5,'alpha',3)
% @endcode
%
for i = 1:length(varargin)/2
   varname = varargin{(i-1)*2+1};
   varvalue = varargin{(i-1)*2+2};
   assignin('caller',varname,varvalue);
end
end

