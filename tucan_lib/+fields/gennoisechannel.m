function [ux uy uz] = gennoisechannel(xux,yux,xuy,yuy,xuz,yuz)
%
% channel with y-walls and forcing in x direction
%
%
%
uxmax = 15.;
noisefactor = 0.15;
devnoise = uxmax*noisefactor;
uxmeanfun = @(y) uxmax*(1-y.^2);
noisefun = @(m,n) devnoise*randn(m,n);

ux = zeros(length(xux), length(yux));
uy = zeros(length(xuy), length(yuy));
uz = zeros(length(xuz), length(zuy));

% ux
for i = 1:length(xux)
   ux(i,:) = uxmeanfun(yux) + noisefun(1,length(yux));
end
uy(:,:) = noisefun(length(xuy),length(yux));
uz(:,:) = noisefun(length(xuz),length(yuz));


return
end
