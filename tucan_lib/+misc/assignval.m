function assignval(name,val)
% @author M.Moriche (Incredibly hard work)
% @date 17-12-2013 by M.Moriche \n
%       Documented
%
% @brief function to assign values from a known workspace.
%
% @warning If someone knows how to assign values of the current
%          workspace, please send an email to mmoriche@ing.uc3m.es 
%          it will be very much appreciated.
%
assignin('caller',name,val);
end
