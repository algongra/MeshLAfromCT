close all
clear all

% generate eulerian mesh for TUCAN
path = '/data2/mmoriche/IBcode/test';
basename1 = 'channelsphere_1';
fr = Frame(basename1, 'path', path)
fr.readfields('items',{'params','buffers'});

% sizes of the new field
nx = 64;
ny = 64;
nz = 64;

% new mesh
[xux1 yux1 zux1 ...
 xuy1 yuy1 zuy1 ...
 xuz1 yuz1 zuz1 ...
 xp1  yp1  zp1  ] = fields.genmesheu(fr.x0,fr.xf...
                                    ,fr.y0,fr.yf...
                                    ,fr.z0,fr.zf,fr.nx,fr.ny,fr.nz);
% new mesh
[xux2 yux2 zux2 ...
 xuy2 yuy2 zuy2 ...
 xuz2 yuz2 zuz2 ...
 xp2  yp2  zp2  ] = fields.genmesheu(fr.x0,fr.xf...
                                    ,fr.y0,fr.yf...
                                    ,fr.z0,fr.zf,nx,ny,nz);

[X1 Y1 Z1] = ndgrid(xux1,yux1,zux1);
[X2 Y2 Z2] = ndgrid(xux2,yux2,zux2);
ux = interpn(X1,Y1,Z1,fr.ux,X2,Y2,Z2);
[X1 Y1 Z1] = ndgrid(xuy1,yuy1,zuy1);
[X2 Y2 Z2] = ndgrid(xuy2,yuy2,zuy2);
uy = interpn(X1,Y1,Z1,fr.uy,X2,Y2,Z2);
[X1 Y1 Z1] = ndgrid(xuz1,yuz1,zuz1);
[X2 Y2 Z2] = ndgrid(xuz2,yuz2,zuz2);
uz = interpn(X1,Y1,Z1,fr.uz,X2,Y2,Z2);
[X1 Y1 Z1] = ndgrid(xp1,yp1,zp1);
[X2 Y2 Z2] = ndgrid(xp2,yp2,zp2);
p = interpn(X1,Y1,Z1,fr.p,X2,Y2,Z2);
[X1 Y1 Z1] = ndgrid(xp1,yp1,zp1);
[X2 Y2 Z2] = ndgrid(xp2,yp2,zp2);
phi = interpn(X1,Y1,Z1,fr.phi,X2,Y2,Z2);


fields.genh5frame(ux,uy,uz,p,phi,'channelsphere_3.h5res');
