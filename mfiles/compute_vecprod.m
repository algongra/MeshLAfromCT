function [area,normal,vecprod] = compute_vecprod(vv,cc)
%
ntri = size(cc,1);
nt   = size(vv,3);
%
area = zeros(ntri,1);
normal = zeros(ntri,3);
vecprod = zeros(ntri,3);
%
v1 = zeros(3,1);
v2 = zeros(3,1);
dum1 = zeros(3,1);
%
for it = 1:ntri
 dum1 = squeeze(vv(cc(it,1),:));
 v1 = squeeze(vv(cc(it,2),:))-dum1;
 v2 = squeeze(vv(cc(it,3),:))-dum1;
 v3 = cross(v1,v2);
 modv = norm(v3);
 area(it,:) = 0.5*modv;
 vecprod(it,:) = v3;
 normal(it,:) = v3/modv;
end
