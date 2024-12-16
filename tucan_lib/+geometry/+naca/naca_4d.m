function [xl xu yl yu]=naca_4d(digits,x,c,varargin)
% @author M.Moriche
%
% @date 13-12-2013 by M.Moriche \n
%       Documented and modify defaults assignment
% 
% @brief Function to obtaion lower and uper x,y coordiantes of
%        a two digits NACA airfoil.
% 
% @details
% 
% Mandatory arguments
%  - digits: four items character string defining the airfoil shape
%  - x: coordinates along the chord of the point to be evaluated.
%  - c: chord length
% 
% Optional arguments
%  - alpha: angle in degrees. (Watchout is -AoA) (default = 0)
%  - controlpoint: [xc;yc] control point (used as a coordinate reference)
%                  (default = [0;0])
%  - zerothickness: boolean to indicate if the coefficients are modified
%                   in such a way that the trailing edge thickness is 
%                   revomed. Mainly used for computational purposes.
%                   (default=false)
% 
% Examples:
% 
% @code
% [xl xu yl yu]=naca4digits('0012',0.1,1)
% [xl xu yl yu]=naca4digits('0012',0.1,1,'controlpoint',[-1;0],'alpha',10)
% [xl xu yl yu]=naca4digits('0012',0.1,1,'zerothickness',true)
% @endcode
% 

alpha = 0;
controlpoint = [0;0];
zerothickness = false;
misc.assigndefaults(varargin{:});


m=str2num(digits(1))/100;
p=str2num(digits(2))/10;
t=str2num(digits(3:4))/100;
%
if zerothickness
   yt=t*c/0.2*(0.2969*(x/c)^0.5 ...
              -0.1260*(x/c)^1   ...
              -0.3516*(x/c)^2   ...
              +0.2843*(x/c)^3   ...
              -0.1036*(x/c)^4);
else
   yt=t*c/0.2*(0.2969*(x/c)^0.5 ...
              -0.1260*(x/c)^1   ...
              -0.3516*(x/c)^2   ...
              +0.2843*(x/c)^3   ...
              -0.1015*(x/c)^4);
end
%
if x < p*c
   yc=m*x/p^2*(2*p-x/c);
   dycdx=m/p^2*(2*p-x/c)+m*x/p^2*(-1/c);
else
   yc=m*(c-x)/(1-p)^2*(1+x/c-2*p);
   dycdx=-m/(1-p)^2*(1+x/c-2*p)+m*(c-x)/(1-p)^2*1/c;
end
%
theta=atan(dycdx);
%
xl=x+yt*sin(theta);
xu=x-yt*sin(theta);
%
yl=yc-yt*cos(theta);
yu=yc+yt*cos(theta);
%
A=[cos(alpha) -sin(alpha);sin(alpha) cos(alpha)];
posl=controlpoint + A*[xl;yl];
posu=controlpoint + A*[xu;yu];
xl=posl(1);yl=posl(2);
xu=posu(1);yu=posu(2);
%
end
