close all
clear all

patt = '/home/mmoriche/channel/retau180.%spxz'

%x0 = 
%xf = 
%y0 = 
%yf = 
%z0 =
%zf = 



% ux
filename = sprintf(patt, 'u');
[time,x,posy,z,wkux]=fields.readfieldxz(filename);

uxmean = zeros(length(posy),1);
for iy = 1:length(posy)
  auxarray = wkux(iy,:,:);
  uxmean(iy) = mean(auxarray(:));
end
figure()
plot(posy,uxmean)
figure()
dudy = diff(uxmean)./diff(posy);
plot(posy(1:end-1),dudy)
