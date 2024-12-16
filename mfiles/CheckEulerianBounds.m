function [x0 y0 z0 xf yf zf nx ny nz] = CheckEulerianBounds(x0,y0,z0,xf,yf,zf,nx,ny,nz,mesh,ib_mesh,dh)

% ACHTUNG!!!
% This function assumes that the last edge is the MV
iMV = mesh.n_edges;

% Check that bounds of the Eulerian grid contain all Lagrangian points
xibmax = max(max(ib_mesh.pts(:,1,:)));
xibmin = min(min(ib_mesh.pts(:,1,:)));
%
yibmax = max(max(ib_mesh.pts(:,2,:)));
yibmin = min(min(ib_mesh.pts(:,2,:)));
%
zibmax = max(max(ib_mesh.pts(:,3,:)));
zibmin = min(min(ib_mesh.pts(:,3,:)));
%
if xibmin < x0 | xibmax > xf |...
   yibmin < y0 | yibmax > yf |...
   zibmin < z0 | zibmax > zf
   wartxt = '\nWarning: Lagrangian points are outside of Eulerian bounds';
   disp(sprintf(wartxt));
   disp(sprintf(['   x0      xibmin     xibmax        xf\n'...
                 '%8.4f   %8.4f   %8.4f   %8.4f'],x0,xibmin,xibmax,xf));
   disp(sprintf(['   y0      yibmin     yibmax        yf\n'...
                 '%8.4f   %8.4f   %8.4f   %8.4f'],y0,yibmin,yibmax,yf));
   disp(sprintf(['   z0      zibmin     zibmax        zf\n'...
                 '%8.4f   %8.4f   %8.4f   %8.4f'],z0,zibmin,zibmax,zf));
   % Update Eulerian bounds values
   x0ori = x0; xfori = xf; y0ori = y0; yfori = yf; z0ori = z0; zfori = zf;
   if xibmin < x0; x0 = floorval(xibmin,2); end;
   if xibmax > xf; xf =  ceilval(xibmax,2); end;
   if yibmin < y0; y0 = floorval(yibmin,2); end;
   if yibmax > yf; yf =  ceilval(yibmax,2); end;
   if zibmin < z0; z0 = floorval(zibmin,2); end;
   if zibmax > zf; zf =  ceilval(zibmax,2); end;
   %
   % Update nx, ny, and nz (implies moving a little bit the bounds)
   disp(sprintf('Eulerian bounds are going to be updated to:'));
   %
   % ACHTUNG!!!
   % We assume that lower/upper limits take always negative/positive values
   if xibmin < x0ori | xibmax > xfori
      nx = (xf-x0)/dh + 1;
      nx0 = ceil(nx*abs(x0)/(xf-x0)); nxf = ceil(nx*abs(xf)/(xf-x0));
      nx = nx0 + nxf; x0 = -nx0*dh; xf = nxf*dh;
      disp(sprintf('x0: %g, xf: %g, nx: %i',x0,xf,nx));
   end
   %
   if yibmin < y0ori | yibmax > yfori
      ny = (yf-y0)/dh + 1;
      ny0 = ceil(ny*abs(y0)/(yf-y0)); nyf = ceil(ny*abs(yf)/(yf-y0));
      ny = ny0 + nyf; y0 = -ny0*dh; yf = nyf*dh;
      disp(sprintf('y0: %g, yf: %g, ny: %i',y0,yf,ny));
   end
   %
   if zibmin < z0ori | zibmax > zfori
      nz = (zf-z0)/dh + 1;
      nz0 = ceil(nz*abs(z0)/(zf-z0)); nzf = ceil(nz*abs(zf)/(zf-z0));
      nz = nz0 + nzf; z0 = -nz0*dh; zf = nzf*dh;
      disp(sprintf('z0: %g, zf: %g, nz: %i',z0,zf,nz));
   end
end


% Check if flow from MV is separated from Eulerian bounds a characteristic
% length in Mitral valve normal direction
%
% Obtain characteristic length of Mitral valve
iMV = mesh.n_edges; ifacMV = find(mesh.face_type==iMV);
LMV = max(sqrt(sum(mesh.areas_mov(ifacMV,:),1)));
iibMV = find(ib_mesh.type==iMV);
%
for ifr = 1:mesh.nt
    % Rename variables
    normalsMV = ib_mesh.normals(iibMV,:,ifr);
    ptsMV = ib_mesh.pts(iibMV,:,ifr);
    %
    normalsMV = normalsMV./VecNorm(normalsMV,2,2);
    auxibpts = ptsMV + LMV*normalsMV;
    %
    xibMVmin(ifr) = min(auxibpts(:,1)); xibMVmax(ifr) = max(auxibpts(:,1));
    yibMVmin(ifr) = min(auxibpts(:,2)); yibMVmax(ifr) = max(auxibpts(:,2));
    zibMVmin(ifr) = min(auxibpts(:,3)); zibMVmax(ifr) = max(auxibpts(:,3));
end
%
if min(xibMVmin) < x0 | max(xibMVmax) > xf |...
   min(yibMVmin) < y0 | max(yibMVmax) > yf |...
   min(zibMVmin) < z0 | max(zibMVmax) > zf
   wartxt = '\nWarning: MV Lagrangian points are too close to Eulerian bounds';
   disp(sprintf(wartxt));
   disp(sprintf(['   x0      xibMVmin   xibMVmax      xf\n'...
                 '%8.4f   %8.4f   %8.4f   %8.4f'],...
                 x0,min(xibMVmin),max(xibMVmax),xf));
   disp(sprintf(['   y0      yibMVmin   yibMVmax      yf\n'...
                 '%8.4f   %8.4f   %8.4f   %8.4f'],...
                 y0,min(yibMVmin),max(yibMVmax),yf));
   disp(sprintf(['   z0      zibMVmin   zibMVmax      zf\n'...
                 '%8.4f   %8.4f   %8.4f   %8.4f'],...
                 z0,min(zibMVmin),max(zibMVmax),zf));
   % Update Eulerian bounds values
   x0ori = x0; xfori = xf; y0ori = y0; yfori = yf; z0ori = z0; zfori = zf;
   if min(xibMVmin) < x0; x0 = floorval(min(xibMVmin),2); end;
   if max(xibMVmax) > xf; xf =  ceilval(max(xibMVmax),2); end;
   if min(yibMVmin) < y0; y0 = floorval(min(yibMVmin),2); end;
   if max(yibMVmax) > yf; yf =  ceilval(max(yibMVmax),2); end;
   if min(zibMVmin) < z0; z0 = floorval(min(zibMVmin),2); end;
   if max(zibMVmax) > zf; zf =  ceilval(max(zibMVmax),2); end;
   %
   % Update nx, ny, and nz (implies moving a little bit the bounds)
   disp(sprintf('Eulerian bounds are going to be updated to:'));
   %
   % ACHTUNG!!!
   % We assume that lower/upper limits take always negative/positive values
   if min(xibMVmin) < x0ori | max(xibMVmax) > xfori
      nx = (xf-x0)/dh + 1;
      nx0 = ceil(nx*abs(x0)/(xf-x0)); nxf = ceil(nx*abs(xf)/(xf-x0));
      nx = nx0 + nxf; x0 = -nx0*dh; xf = nxf*dh;
      disp(sprintf('x0: %g, xf: %g, nx: %i',x0,xf,nx));
   end
   %
   if min(yibMVmin) < y0ori | max(yibMVmax) > yfori
      ny = (yf-y0)/dh + 1;
      ny0 = ceil(ny*abs(y0)/(yf-y0)); nyf = ceil(ny*abs(yf)/(yf-y0));
      ny = ny0 + nyf; y0 = -ny0*dh; yf = nyf*dh;
      disp(sprintf('y0: %g, yf: %g, ny: %i',y0,yf,ny));
   end
   %
   if min(zibMVmin) < z0ori | max(zibMVmax) > zfori
      nz = (zf-z0)/dh + 1;
      nz0 = ceil(nz*abs(z0)/(zf-z0)); nzf = ceil(nz*abs(zf)/(zf-z0));
      nz = nz0 + nzf; z0 = -nz0*dh; zf = nzf*dh;
      disp(sprintf('z0: %g, zf: %g, nz: %i',z0,zf,nz));
   end
end

return

end

