function gen3D_z_interp(self, kz) 

nz = self.nz;
z0 = self.z0;
zf = self.zf;

ux = zeros([size(self.frFLO.ux_re)    nz]);
uy = zeros([size(self.frFLO.uy_re)    nz]);
uz = zeros([size(self.frFLO.uz_re)  nz-1]);
p  = zeros([size(self.frFLO.p_re)     nz]);

[XXux YYux]   = ndgrid(self.fr2D.xux, self.fr2D.yux);
[XXuy YYuy]   = ndgrid(self.fr2D.xuy, self.fr2D.yuy);
[XXp YYp]     = ndgrid(self.fr2D.xp , self.fr2D.yp);
[XXfux YYfux] = ndgrid(self.frFLO.xux, self.frFLO.yux);
[XXfuy YYfuy] = ndgrid(self.frFLO.xuy, self.frFLO.yuy);
[XXfp YYfp]   = ndgrid(self.frFLO.xp , self.frFLO.yp);

zux = linspace(z0, zf, nz);
zuy = linspace(z0, zf, nz);
zp  = linspace(z0, zf, nz);
zuz = linspace(z0, zf, nz);
zuz = zuz(1:end-1)+self.fr2D.dx/2;

b1 = interpn(XXux, YYux, self.fr2D.ux, XXfux, YYfux);
b2 = interpn(XXuy, YYuy, self.fr2D.uy, XXfuy, YYfuy);
b3 = interpn(XXp , YYp , self.fr2D.p , XXfp , YYfp);

for i0 = 1:nz

   zi = zux(i0);

   self.frFLO.C2R_z(kz, zi, 'varlist', {'ux','uy','p'});

   c = self.frFLO.ux;
   a = b1+c;
   ux(:,:,i0) =a;

   c = self.frFLO.uy;
   a = b2+c;
   uy(:,:,i0) =a;

   c = self.frFLO.p;
   a = b3+c;
   p(:,:,i0) =a;
end

for i0 = 1:nz-1

   zi = zuz(i0);

   self.frFLO.C2R_z(kz, zi, 'varlist', {'uz'});

   c = self.frFLO.uz;
   uz(:,:,i0) =c;

end

self.ux = ux;
self.uy = uy;
self.uz = uz;
self.p  = p;

self.zux = zux;
self.zuy = zuy;
self.zuz = zuz;
self.zp  = zp;

self.xux = self.frFLO.xux;
self.xuy = self.frFLO.xuy;
self.xuz = self.frFLO.xuz;
self.xp  = self.frFLO.xp;

self.yux = self.frFLO.yux;
self.yuy = self.frFLO.yuy;
self.yuz = self.frFLO.yuz;
self.yp  = self.frFLO.yp;

self.sizebyvar.ux.count = size(ux);
self.sizebyvar.uy.count = size(uy);
self.sizebyvar.uz.count = size(uz);
self.sizebyvar.p.count  = size(p);

return
end
