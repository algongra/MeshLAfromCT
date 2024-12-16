function vorticity(self,dim)
% @author M.Moriche
% @date Modified 25-05-2013 by M.Moriche
% @date Modified 06-05-2016 by A.Gonzalo
% @date 26-07-2016 by M.Moriche \n
%       bug in vorticity(2), gave -wy as output
% @date 26-07-2016 by M.Moriche \n
%       bug in vorticity(2), gave -wy as output
%
% @brief  Calculates vorticity
%    
% @details
%
%
% Calculates vorticity at xux, yuy points.
% Is saved in the property "w"
%
% @verbatim
%
% o  -  o  -  o  -  o  -  o  -  o      o  pressure points
%                                      - ux points
% |  *  |  *  |  *  |  *  |  *  |      | uy points
%                                
% o  -  o  -  o  -  o  -  o  -  o      
%                                      * vorticity points
% |  *  |  *  |  *  |  *  |  *  |
%                                
% o  -  o  -  o  -  o  -  o  -  o 
%                                
% |  *  |  *  |  *  |  *  |  *  |
%                                
% o  -  o  -  o  -  o  -  o  -  o 
%
% @endverbatim
%
%
% Usage example:
%   
% @code
% fr = SUFF(2);
% fr.vorticity(3);
% imagesc(fr.wz')
% @endcode
 
coords = {'x','y','z'};

dir1ref = [2 3 1];
dir2ref = [3 1 2];
dir1 = dir1ref(dim);
dir2 = dir2ref(dim);

nm  = coords{dim};
nm1 = coords{dir1};
nm2 = coords{dir2};

fld1 = getfield(self, ['u' nm1]);
fld2 = getfield(self, ['u' nm2]);

c1 = getfield(self, [nm1 'u' nm2]);
c2 = getfield(self, [nm2 'u' nm1]);

dc1 = diff(c1);
dc2 = diff(c2);

% needs to reshape and give dim3 to avoid merging of dims by repmat 
% for example, if dc1 = dz, needs size(dz) = [1 1 nz]
tore = [1 1 1];
tore(dir1) = length(dc1);
dc1 = reshape(dc1,tore);

tore = [1 1 1];
tore(dir2) = length(dc2);
dc2 = reshape(dc2,tore);

i1 = size(diff(fld2,1,dir1));
rep1 = i1;
rep1(dir1) = 1;
h1 = repmat(dc1,rep1);

i2 = size(diff(fld1,1,dir2));
rep2 = i2;
rep2(dir2) = 1;
h2 = repmat(dc2,rep2);

iend = min(i1(1),i2(1));
jend = min(i1(2),i2(2));
kend = min(i1(3),i2(3));

d1 = diff(fld2,1,dir1)./h1;
d2 = diff(fld1,1,dir2)./h2;

self.(['w' nm]) =  d1(1:iend, 1:jend, 1:kend) - d2(1:iend, 1:jend, 1:kend);

end
