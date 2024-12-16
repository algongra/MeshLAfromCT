function [xyz da vol] = spiral(D,dx)

A = 4*pi*(D/2)^2;  % area of the sphere carcase
N = round(A/dx^2);

k=1:N;
h = -1+ 2*(k-1)/(N-1);
theta = acos(h);

k=1;
nextphi = @(oldphi,hk) mod(oldphi + 3.6/(sqrt(N)*sqrt(1-hk^2)),2*pi);

phi(1) = 0;
for k=2:N-1
   phi(k) = nextphi(phi(k-1), h(k));
end
phi(N) = 2*pi;
r = ones(size(phi));

[x y z] = sph2cart(theta, phi, r);

xyz = zeros(length(x),3);
xyz(:,1) = x;
xyz(:,2) = y;
xyz(:,3) = z;

da = A/N;
vol = da*dx;

return
end
