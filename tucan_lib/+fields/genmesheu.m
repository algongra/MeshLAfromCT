function [xux yux zux ...
          xuy yuy zuy ...
          xuz yuz zuz ...
          xp  yp  zp  ] = genmesheu(x0,xf,y0,yf,z0,zf,nx,ny,nz,varargin)

perx = 0;
pery = 0;
perz = 0;
ghost = false;
misc.assigndefaults(varargin{:});

xux = fields.genarray(x0,xf,nx+perx,'stag',1,'ghost',ghost);
yux = fields.genarray(y0,yf,ny,'ghost',ghost);
zux = fields.genarray(z0,zf,nz,'ghost',ghost);

xuy = fields.genarray(x0,xf,nx,'ghost',ghost);
yuy = fields.genarray(y0,yf,ny+pery,'stag',1,'ghost',ghost);
zuy = fields.genarray(z0,zf,nz,'ghost',ghost);

xuz = fields.genarray(x0,xf,nx,'ghost',ghost);
yuz = fields.genarray(y0,yf,ny,'ghost',ghost);
zuz = fields.genarray(z0,zf,nz+perz,'stag',1,'ghost',ghost);

xp  = fields.genarray(x0,xf,nx,'ghost',ghost);
yp  = fields.genarray(y0,yf,ny,'ghost',ghost);
zp  = fields.genarray(z0,zf,nz,'ghost',ghost);

return
end
