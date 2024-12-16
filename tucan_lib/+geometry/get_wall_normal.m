function [unx,uny,utx,uty] =  get_wall_normal(x,y,varargin)
   %
   %
   % Mandatory arguments:
   %  - x
   %  - x
   %
   % Optional arguments:
   %  - closed: boolean to indicate if the curve is closed or open
   %
   closed=true;
   misc.assigndefaults(varargin{:});
   N = length(x);

   vx = zeros(size(x));
   vy = zeros(size(y));
   vx(2:end-1) = x(3:end) - x(1:end-2);
   vy(2:end-1) = y(3:end) - y(1:end-2);
   if closed
      vx(1)  = x(2) - x(end);
      vx(end) = x(1) - x(end-1);
      vy(1)  = y(2) - y(end);
      vy(end) = y(1) - y(end-1);
   else
      vx(1)  = x(2)  - x(1);
      vx(end) = x(end) - x(end-1);
      vy(1)  = y(2)  - y(1);
      vy(end) = y(end) - y(end-1);
   end

   vmod = sqrt( (vx.^2 + vy.^2) );
   utx = vx./vmod;
   uty = vy./vmod;

   unx = zeros(size(x));
   uny = zeros(size(y));
   for i = 1:N
       A =  [utx(i)  uty(i);...
             uty(i) -utx(i)];
       b =   [0 ; 1];
       sol = A\b;
       unx(i) = sol(1);
       uny(i) = sol(2);
   end
end
