function [xux yux zux ...
          xuy yuy zuy ...
          xuz yuz zuz ...
          xp  yp  zp  ] = mesheu(x0,xf,y0,yf,z0,zf,nx,ny,nz,varargin)

% ARGUMENTS USED TO GET THIS MESH:
% mandatory arguments:
% x0 = -1; xf = 1; nx = 128;
%
% y0 = -1; yf = 1; ny = 128;
%
% z0 = -1; zf = 1; nz = 128;
%
% varargin arguments:
% 'perx',1,'pery',1,'perz',1,'ghost',true,'write',true
%
% addpath('~/source/IBcode/lib/mfiles'); 

perx = 0;
pery = 0;
perz = 0;
bcfixx = false;
bcfixy = false;
bcfixz = false;
ghost = false;
write = false;
misc.assigndefaults(varargin{:});

%------------------------------------------------------------------------------% 
%
%
%                  ## ### ###     ##  # #     # #  ## ### ##  
%                 #   #    #      # # # #     # # #   #   # # 
%                  #  ##   #      ##   #      # #  #  ##  ##  
%                   # #    #      # #  #      # #   # #   # # 
%                 ##  ###  #      ##   #      ### ##  ### # # 
%
%
%            -- CHOOSE THE FUNCTION TO CREATE x DIRECTION MESH --
%
% Set parameters needed 
dir = 1;
%
[xux xuy xuz xp dxr nx irefbeg] = fields.genuniform(x0,xf,dir,nx,...
                                                    'ghost',ghost,'per',perx);
%
%            -- CHOOSE THE FUNCTION TO CREATE y DIRECTION MESH --
%
% Set parameters needed 
dir = 2;
%
[yux yuy yuz yp dyr ny jrefbeg] = fields.genuniform(y0,yf,dir,ny,...
                                                    'ghost',ghost,'per',pery);
%
%            -- CHOOSE THE FUNCTION TO CREATE z DIRECTION MESH --
%
% Set parameters needed 
dir = 3;
%
[zux zuy zuz zp dzr nz krefbeg] = fields.genuniform(z0,zf,dir,nz,...
                                                    'ghost',ghost,'per',perz);
%------------------------------------------------------------------------------%

fprintf('\nmesh generated')

if write
   % Set number of points in each direction as integers (nx, ny, nz)
   nx = int32(nx); ny = int32(ny); nz = int32(nz);
   % Set beggining of refined zone indices arrays as integers
   irefbeg = int32(irefbeg); jrefbeg = int32(jrefbeg); krefbeg = int32(krefbeg);

   fext = '.h5mesh';
   funpth = mfilename('fullpath'); funnm = mfilename;
   filenm = strcat(funnm,fext); fullnm = strcat(funpth,fext);

   if ghost
      % Change the orientation (c arrays instead of fortran arrays)
      xux0 = xux'; yux0 = yux'; zux0 = zux'; xp0 = xp';
      xuy0 = xuy'; yuy0 = yuy'; zuy0 = zuy'; yp0 = yp';
      xuz0 = xuz'; yuz0 = yuz'; zuz0 = zuz'; zp0 = zp';
   else
      fprintf('\nError writting mesh into HDF5 file in %s function\n',funnm)
      fprintf('IBcode need meshes with ghost points\n')
      fprintf('Call %s adding argument ''ghost'',true\n\n',funnm)
      quit
   end

   meshlist = {'nx','xux0','yux0','zux0','xp0','dxr','irefbeg',...
               'ny','xuy0','yuy0','zuy0','yp0','dyr','jrefbeg',...
               'nz','xuz0','yuz0','zuz0','zp0','dzr','krefbeg'};

   io.writeh5file(meshlist,fullnm);

   fprintf('\n%s file generated\n\n',filenm)
end

return
end
