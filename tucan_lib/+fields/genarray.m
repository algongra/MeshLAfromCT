function x = genarray(x0,xf,nx,varargin)
% @author M.Moriche
% @date 10-10-2014 by M.Moriche \

stag = 0; % stag  can be 0 or 1
per = 0;
ghost = false;
misc.assigndefaults(varargin{:});

dx = (xf-x0)/double(nx);

n = double(nx) - stag + per;
if ~ghost
   x = x0 + dx*(1:n) - (1-stag)*0.5*dx;
else
   x = x0 + dx*(0:n+1) - (1-stag)*0.5*dx;
end
%x = linspace(x0,xf-(1-stag)*0.5*dx,n);

return
end
