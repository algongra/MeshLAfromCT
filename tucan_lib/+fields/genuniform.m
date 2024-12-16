function [cux cuy cuz cp dr n indxrb] = genuniform(c0,cf,dir,n,varargin)

% @author A.Gonzalo
% 
% @brief function to build one direction of a 3D mesh with the same
%        resolution in the whole mesh
% 
% @date  20-01-2015 by A.Gonzalo\n
%        Created
% @date  20-01-2015 by A.Gonzalo\n
%        Documented
% @date  09-04-2015 by A.Gonzalo\n
%        Added indxrb
% 
% @details
%
% + MANDATORY ARGUMENTS:
%  --------------------
%
%  - c0: initial physical point of the mesh [ghost point mesh]
%  - cf: final physical point of the mesh [ghost point mesh]
%  - dir: direction (x=1, y=2, z=3)
%  - n: number of point in the mesh [solver points of pressure mesh]
%
% + OUTPUT ARGUMENTS:
%  -----------------
%
%  - cux: x velocity points [face/body centered mesh]
%  - cuy: y velocity points [face/body centered mesh]
%  - cuz: z velocity points [face/body centered mesh]
%  - cp: pressure points [body centered mesh]
%  - dc: distance between points (dxr, dyr, dzr)
%  - n: number of point in the mesh [solver points of pressure mesh]
%  - indxrb: beggining of refined zone index
%
%
% @code
% dir=1;
% [xux xuy xuz xp dxr irefbeg] = fields.genuniform(c0,cf,dir,n,'dr',0.01);
%
% dir=2;
% [yux yuy yuz yp dyr jrefbeg] = fields.genuniform(c0,cf,dir,n,'per',1);
%
% dir=3;
% [zux0 zuy0 zuz0 zp0 dzr krefbeg] = fields.genuniform(c0,cf,dir,n,...
%                                              'ghost',true,'per',1);
% @endcode

dr = 0;
per = 0;
ghost = false;
misc.assigndefaults(varargin{:});

if dr ~= 0
   n = round((cf-c0)/dr);
   if per == 1
      nold = n;
      naux = nextpow2(n);
      nhi = 2^naux;
      nlo = 2^(naux-1);
      if (nhi-n) <= (n-nlo)
         n = nhi;
      else
         n = nlo;
      end
      % check number of points for periodic mesh
      if nold ~= n
         drnlo = (cf-c0)/nlo;
         drnhi = (cf-c0)/nhi;
         lengnlo = dr*nlo;
         lengnhi = dr*nhi;
         funcnm = mfilename;
         disp(' ')
         disp(' ')
         fprintf('ERROR SETTING MESH PARAMETERS IN %s FUNCTION',funcnm)
         disp(' ')
         disp(' ')
         aggfrmt = ['With this value of dr and those bounds values (c0 and '...
                    'cf) HYPRE cannot solve a PERIODIC mesh'];
         disp(aggfrmt)
         aggfrmt = ['You can use the same bounds using dr = %18.15f to get '...
                    'a mesh with n = %i or dr = %18.15f to get a mesh with '...
                    'n = %i\n'];
         fprintf(aggfrmt,drnlo,nlo,drnhi,nhi);
         aggfrmt = ['or use the same dr using (cf-c0) = %18.15f to get a '...
                    'mesh with n = %i or (cf-c0) = %18.15f to get a mesh '...
                    'with n = %i\n'];
         fprintf(aggfrmt,lengnlo,nlo,lengnhi,nhi);
         disp(' ')
         disp(' ')
         % give values to the output arguments to avoid the MATLAB error display
         cux=0; cuy=0; cuz=0; cp=0; dr=0; n=0;
         quit
      end
   end
end

dr = (cf-c0)/n;
indxrb = zeros(3,1);

if dir == 1
   cux = fields.genarray(c0,cf,n,'stag',1,'ghost',ghost,'per',per);
   cuy = fields.genarray(c0,cf,n,'ghost',ghost);
   cuz = fields.genarray(c0,cf,n,'ghost',ghost);
   cp  = fields.genarray(c0,cf,n,'ghost',ghost);
elseif dir == 2
   cux = fields.genarray(c0,cf,n,'ghost',ghost);
   cuy = fields.genarray(c0,cf,n,'stag',1,'ghost',ghost,'per',per);
   cuz = fields.genarray(c0,cf,n,'ghost',ghost);
   cp  = fields.genarray(c0,cf,n,'ghost',ghost);
elseif dir == 3
   cux = fields.genarray(c0,cf,n,'ghost',ghost);
   cuy = fields.genarray(c0,cf,n,'ghost',ghost);
   cuz = fields.genarray(c0,cf,n,'stag',1,'ghost',ghost,'per',per);
   cp  = fields.genarray(c0,cf,n,'ghost',ghost);
end

return
end
