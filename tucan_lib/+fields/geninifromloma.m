close all
clear all

patt = '/home/mmoriche/channel/cont99.%spxz'

%x0 = 
%xf = 
%y0 = 
%yf = 
%z0 =
%zf = 



% ux
filename = sprintf(patt, 'u');
[time,x,posy,z,wkux]=fields.readfieldxz(filename);

%uxmean = zeros(length(posy),1);
%for iy = 1:length(posy)
   %auxarray = wkux(iy,:,:);
   %uxmean(iy) = mean(auxarray(:));
%end
%figure()
%plot(posy,uxmean)
%figure()
%dudy = diff(uxmean)./diff(posy);
%plot(posy(1:end-1),dudy)


% uy
filename = sprintf(patt, 'v');
[time,x,posy,z,wkuy]=fields.readfieldxz(filename);

% uz
filename = sprintf(patt, 'w');
[time,x,posy,z,wkuz]=fields.readfieldxz(filename);

% scale fields with 1/u_tau approx = 30
wkux = 15.71014.*wkux;
wkuy = 15.71014.*wkuy;
wkuz = 15.71014.*wkuz;


% generate eulerian mesh for TUCAN

x0 = x(1);
xf = x(end);
y0 = posy(1);
yf = posy(end);
z0 = z(1);
zf = z(end);

nx = 256;
ny = 360;
nz = 180;

[xux yux zux ...
 xuy yuy zuy ...
 xuz yuz zuz ...
 xp  yp  zp  ] = fields.genmesheu(x0,xf,y0,yf,z0,zf,nx,ny,nz);

[X  Y  Z ] = meshgrid(x,posy,z);
[X2 Y2 Z2] = meshgrid(xux,yux,zux);
ux = interp3(X,Y,Z,wkux,X2,Y2,Z2);
[X2 Y2 Z2] = meshgrid(xuy,yuy,zuy);
uy = interp3(X,Y,Z,wkuy,X2,Y2,Z2);
[X2 Y2 Z2] = meshgrid(xuz,yuz,zuz);
uz = interp3(X,Y,Z,wkuz,X2,Y2,Z2);


ux2 = zeros(length(xux), length(yux), length(zux));
uy2 = zeros(length(xuy), length(yuy), length(zuy));
uz2 = zeros(length(xuz), length(yuz), length(zuz));

% flip dimensions
for i = 1:length(xux)
   for j = 1:length(yux)
      for k = 1:length(zux)
         ux2(i,j,k) = ux(j,i,k);
      end
   end
end

for i = 1:length(xuy)
   for j = 1:length(yuy)
      for k = 1:length(zuy)
         uy2(i,j,k) = uy(j,i,k);
      end
   end
end

for i = 1:length(xuz)
   for j = 1:length(yuz)
      for k = 1:length(zuz)
         uz2(i,j,k) = uz(j,i,k);
      end
   end
end

fields.genh5frame(ux2,uy2,uz2,'channeltoni.h5frame');
