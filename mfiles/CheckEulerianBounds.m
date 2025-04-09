function [x0 y0 z0 xf yf zf nx ny nz] = CheckEulerianBounds(x0,y0,z0,xf,yf,zf,nx,ny,nz,mesh,ib_mesh,dh)

% nvalid expected for different resolutions
%
% lowres  (first  column)
% nomres  (second column)
% highres (third  column)
nvalid = [144 150 160 162 168 180 192 196 200 210 216 224 240 250 252 ...
          256 270 280 288 294 300 320 324 350 360 378 ...
	  384 392 400 420 432 448 450 480 486 490 500 504 512];

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


% Check if length of interior points of coordinate arrays (nx, ny, and nz) is 
% compatible with spectral solver:
%  - ftt requires coordinate arrays to have a length with at least one factor
%    equal to two.
%  - fft works better when coordinates arrays have a length that is power of 2
%    (default is nx=ny=nz=[144|256} for [low|nom] resolution).
%    It accepts 3, 5, or 7 as factors before stop working, but there is a
%    penalty in speed up.
% Thus, when bounds have grown with respect to the default values, they must be
% expand to fulfill these requirements (and length has to be even to ensure 2 is
% alwawys included)
%
% Find values greater than or equal to nx, ny, and nz in nvalid
nx_gr = nvalid(nvalid > nx);
ny_gr = nvalid(nvalid > ny);
nz_gr = nvalid(nvalid > nz);
%
% nx
if ~ismember(nx,nvalid)
   if isempty(nx_gr)
       error('No value in nvalid is greater than or equal to the nx');
   else
       disp(sprintf('\nWarning: nx factors are > 7: Spectral solver issue'));
       disp(sprintf('Eulerian bounds in x direction are going to be updated to:'));
       nx_add = min(nx_gr) - nx;
       if abs(x0) > abs(xf)
          % Update x0
          nx_add_left = ceil(nx_add/2); x0 = x0 - nx_add_left*dh;
          % Update xf
          nx_add_right = floor(nx_add/2); xf = xf + nx_add_right*dh;
       else
          % Update x0
          nx_add_left = floor(nx_add/2); x0 = x0 - nx_add_left*dh;
          % Update xf
          nx_add_right = ceil(nx_add/2); xf = xf + nx_add_right*dh;
       end
       nx = nx + nx_add;
       disp(sprintf('x0: %g, xf: %g, nx: %i',x0,xf,nx));
   end
end
% ny
if ~ismember(ny,nvalid)
   if isempty(ny_gr)
       error('No value in nvalid is greater than or equal to the ny');
   else
       disp(sprintf('\nWarning: ny factors are > 7: Spectral solver issue'));
       disp(sprintf('Eulerian bounds in y direction are going to be updated to:'));
       ny_add = min(ny_gr) - ny;
       if abs(y0) > abs(yf)
          % Update y0
          ny_add_left = ceil(ny_add/2); y0 = y0 - ny_add_left*dh;
          % Update yf
          ny_add_right = floor(ny_add/2); yf = yf + ny_add_right*dh;
       else
          % Update y0
          ny_add_left = floor(ny_add/2); y0 = y0 - ny_add_left*dh;
          % Update yf
          ny_add_right = ceil(ny_add/2); yf = yf + ny_add_right*dh;
       end
       ny = ny + ny_add;
       disp(sprintf('y0: %g, yf: %g, ny: %i',y0,yf,ny));
   end
end
% nz
if ~ismember(nz,nvalid)
   if isempty(nz_gr)
       error('No value in nvalid is greater than or equal to the nz');
   else
       disp(sprintf('\nWarning: nz factors are > 7: Spectral solver issue'));
       disp(sprintf('Eulerian bounds in y direction are going to be updated to:'));
       nz_add = min(nz_gr) - nz;
       if abs(z0) > abs(zf)
          % Update z0
          nz_add_left = ceil(nz_add/2); z0 = z0 - nz_add_left*dh;
          % Update yf
          nz_add_right = floor(nz_add/2); zf = zf + nz_add_right*dh;
       else
          % Update z0
          nz_add_left = floor(nz_add/2); z0 = z0 - nz_add_left*dh;
          % Update yf
          nz_add_right = ceil(nz_add/2); zf = zf + nz_add_right*dh;
       end
       nz = nz + nz_add;
       disp(sprintf('z0: %g, zf: %g, nz: %i',z0,zf,nz));
   end
end


return

end

