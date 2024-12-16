function gen3D_z(self, kz) 

nz = self.nz;
z0 = self.z0;
zf = self.zf;

ux = zeros([size(self.frFLO.ux_re)    nz]);
uy = zeros([size(self.frFLO.uy_re)    nz]);
uz = zeros([size(self.frFLO.uz_re)  nz-1]);
p  = zeros([size(self.frFLO.p_re)     nz]);

zux = linspace(z0, zf, nz);
zuy = linspace(z0, zf, nz);
zp  = linspace(z0, zf, nz);
zuz = linspace(z0, zf, nz);
zuz = zuz(1:end-1)+self.fr2D.dx/2;

for i0 = 1:nz

   zi = zux(i0);

   self.frFLO.C2R_z(kz, zi, {'ux','uy','p'});

   b = self.fr2D.ux;
   c = self.frFLO.ux;
   a = b+c;
   ux(:,:,i0) =a;

   b = self.fr2D.uy;
   c = self.frFLO.uy;
   a = b+c;
   uy(:,:,i0) =a;

   b = self.fr2D.p;
   c = self.frFLO.p;
   a = b+c;
   p(:,:,i0) =a;
end

for i0 = 1:nz-1

   zi = zux(i0);

   self.frFLO.C2R_z(self.frFLO.kz, zi, {'uz'});

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
