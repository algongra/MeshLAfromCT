function [uproj vproj] = project_quantity(u,v,un,vn)
   mod = sqrt(un.^2 + vn.^2);
   un = un./mod;
   vn = vn./mod;
   dotprod = un.*u + vn.*v;
   uproj = dotprod.*un;
   vproj = dotprod.*vn;
end
