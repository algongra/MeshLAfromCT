function [xy varargout] = circle(D, dx)
% @author M.Moriche
% @brief Function to generate points laying on a circumference.
%
% @details
%
% Generates points laying on a circumference. The arclength between points
%  is the closest value possible to dx, having an associated volume
%  of the marker $vol similar to dx^2.
%
% Reference followed Uhlmann2005 (Appendix A.)
%
% MINIMUM OUTPUT
%  - xy: matrix containing the coordinates of the Lagrangian points
%        whose shape is (n,2)
%
% OPTIONAL OUTPUT
%  - ds: arc length between points.
%  - vol: associated marker volume.
%
% MANDATORY ARGUMENTS
%  - D: circle diameter.
%  - dx: reference mesh width.
%
% NO OPTIONAL ARGUMENTS
%
% EXAMPLES:
% @code
% xy = circle(1., 0.01);
% [xy ds vol] = circle(1., 0.01);
% @endcode
%
R=D/2;
n=round(pi*D/dx);
xy = zeros(n,2);

for i = 1:n
    angle = 2*pi/n*(i-1);
    xy(i,1) =  R*cos(angle);
    xy(i,2) =  R*sin(angle);
end

ds = 2*pi*R/n;
vol=ds*dx;

if nargout == 2
   varargout{1} = ds;
elseif nargout == 3
   varargout{1} = ds;
   varargout{2} = vol;
end

return
end
