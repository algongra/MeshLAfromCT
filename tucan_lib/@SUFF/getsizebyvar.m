function sizebyvar = getsizebyvar(self)

slab = self.slab;
ndim = self.ndim;
coords = {'x','y','z'};

% ij(ivar, idim, ipos)
% ipos=1 beginning
% ipos=2 ending
ij = zeros(ndim+1,ndim,2); 
ivar = 1; % set for p

for ipos = 1:2
   for idim = 1:ndim
      xxname = [coords{idim} 'p'];
      %xx = self.(xxname);
      %xx = io.readh5dset(fnm, xxname);
      xx = self.(xxname);
      xi = slab(idim, ipos);
      if isnan(xi)
         if ipos == 1
            ia = 1;
         elseif ipos == 2
            ia = length(xx);
         end
      else
         ia = find(abs(xx - xi) == min(abs(xx - xi)));
      end
      if (ipos == 1)
         ia = ia(1);
      else
         ia = ia(end);
      end
      ij(1,idim,ipos) = ia;
   end
end

for idim = 1:ndim
   ij(1,idim,2) = max(ij(1,idim,2),ij(1,idim,1)+1);
end

staggarray  = zeros(ndim+1,ndim);
for idim = 1:ndim
   stagarray(idim+1,idim) = 1; % only velocities in its own direction
end

for ivar = 2:ndim+1
   for idim = 1:ndim
      st = stagarray(ivar,idim);
      ij(ivar, idim, 1 ) = ij(1,idim,1)+0;
      ij(ivar, idim, 2 ) = ij(1,idim,2)-st;
   end
end

% velocities
for ivar = 1:ndim
   varnm = ['u' coords{ivar}];
   sizebyvar.(varnm).offset = ij(ivar+1,:,1)-1;
   sizebyvar.(varnm).stride =   ones(1,ndim);
   sizebyvar.(varnm).count =  ij(ivar+1,:,2)-ij(ivar+1,:,1)+1;;
   sizebyvar.(varnm).block =    ones(1,ndim);
end
% pressure
ivar = 1;
varnm = 'p';
sizebyvar.(varnm).offset = ij(ivar,:,1)-1;
sizebyvar.(varnm).stride =   ones(1,ndim);
sizebyvar.(varnm).count =  ij(ivar,:,2)-ij(ivar,:,1)+1;;
sizebyvar.(varnm).block =    ones(1,ndim);

% eulerian position arrays for velocities
for ivar = 1:ndim
   for idim = 1:ndim
      varnm = [coords{idim} 'u' coords{ivar}];
      sizebyvar.(varnm).offset = [ij(ivar+1,idim,1)-1];
      sizebyvar.(varnm).stride =                 [1];
      sizebyvar.(varnm).count =  [ij(ivar+1,idim,2)-ij(ivar+1,idim,1)+1];
      sizebyvar.(varnm).block =                  [1];
   end
end
for idim = 1:ndim
   varnm = [coords{idim} 'p'];
   sizebyvar.(varnm).offset = [ij(1,idim,1)-1];
   sizebyvar.(varnm).stride =                 [1];
   sizebyvar.(varnm).count =  [ij(1,idim,2)-ij(1,idim,1)+1];
   sizebyvar.(varnm).block =                  [1];
end

self.sizebyvar = sizebyvar;

return
end
