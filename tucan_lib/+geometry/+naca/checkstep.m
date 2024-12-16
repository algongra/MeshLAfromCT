function  flag = checkstep(xnow, xnext_1, xnext_2,dsmin, minangle);
% @author M.Moriche
%
% @date 13-12-2013 by M.Moriche \n
%       Documented.
%
% @brief Checks that the vectors xnow-->xnext_1 
%        and xnow-->xnext_2 fulfill certain requierements
%
% @details
% 
% Checks that the angle between vectors is lower thatn minangle
% and that the norm of each of them is lower than dsmin
%
 
v1 = xnext_1 - xnow;
v2 = xnext_2 - xnow;

ang = acos(dot(v1,v2)/(norm(v1)*norm(v2)));
ds = max(norm(v1),norm(v2));

flag = ang < minangle && ds < dsmin;


return 
end
