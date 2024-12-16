function [xf varargout] = circlefull(D, dx)
% @author M.Moriche
% @brief Function to generate points laying on a circle.
%
% @details
%
% Generates points laying on a circle. The arclength between points
%  is the closest value possible to dx, having an associated volume
%  of the marker $vol similar to dx^2.
%
% Reference followed Uhlmann2005 (Appendix A.)
%
% MINIMUM OUTPUT
%  - xf: matrix containing the coordinates of the Lagrangian points
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
% xf = circlefull(1., 0.01);
% [xf ds vol] = circlefull(1., 0.01);
% @endcode
%

R=D/2;
nr=round(R/dx);

nmax=10000;
xf = zeros(nmax,2);

k = 1;
for i = 1:nr
    ri = R-(i-1)*dx;
    perimeter = 2*ri*pi;
    ntheta = round(floor(perimeter/dx));
    for j = 1:ntheta
       angle = 2*pi/ntheta*(j-1);
       xf(k,1) =  ri*cos(angle);
       xf(k,2) =  ri*sin(angle);
       k=k+1;
    end
end
n = k-1;
xf = xf(1:n,:);

perimeter = 2*pi*R;
ntheta = round(floor(perimeter/dx));
ds = 2*pi*R/ntheta;
vol=ds*dx;

if nargout == 2
   varargout{1} = ds;
elseif nargout == 3
   varargout{1} = ds;
   varargout{2} = vol;
end

return
end
