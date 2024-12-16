classdef Frame3D < dynamicprops & Frame
%
% @date 03-07-2014 by O. Flores \n 
%       Added secondinvariant method
%
%
%
% @details       
% methods:       
%  -secondinvariant
%

properties
end
properties(Access='protected')
end
methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function self = Frame3D(basenm,varargin)

   slab = [NaN,NaN,NaN,NaN,NaN,NaN];
   misc.assigndefaults(varargin{:});

   ndim = 3;
   self@Frame(basenm,ndim,varargin{:});

end

function lambda2plot(self,varargin)

   overrides = true;
   misc.assigndefaults(varargin{:});

   if overrides
      self.getlambda();
   elseif isempty(find(strcmp(fields(self), 'l2')))
      self.getlambda();
   end

   self.readfields('items',{'mesheu','meshlag'},'overrides',overrides);

   [Y X Z] = meshgrid(self.yp(2:end-1), self.xp(2:end-1), self.zp(2:end-1));

   l2thres = 0;
   h1 = myplot.isosurfaceplot(X,Y,Z,self.l2,X,l2thres);
   hold on
   plot3(self.yib, self.xib, self.zib, 'k.')
   colormap(jet);
   axis equal
   view(3)
   camroll(-110)
   
return
end

function secondinvariantplot(self, varargin) 

   overrides = true;
   misc.assigndefaults(varargin{:});

   if overrides
      self.secondinvariant();
   elseif isempty(find(strcmp(fields(self), 'Q')))
      self.secondinvariant();
   end

   self.readfields('items',{'mesheu','meshlag'},'overrides',overrides);

   [Y X Z] = meshgrid(self.yp(2:end-1), self.xp(2:end-1), self.zp(2:end-1));


   Qthres = mean(self.Q(:)) + 0.5*std(self.Q(:));
   [f v] = isosurface(Y,X,Z,self.Q, Qthres);

   h1 = patch('Faces',f,'Vertices',v);
   isonormals(Y,X,Z,self.Q,h1);
   set(h1,'facecolor','interp','edgecolor','none');
   set(gca, 'XDir', 'reverse');

   hold on
   plot3(self.yib, self.xib, self.zib, 'k.')

   isocolors(Y,X,Z,X,h1);
   colormap(jet);
   axis equal
   view(3)
   camroll(-110)

return
end

function secondinvariant(self, varargin) 
   % @author O. FLores
   % @date 03-07-2014 by O. Flores \n 
   %       Created
   % @date 13-12-2014 by M. Moriche \n 
   %       Modified. Q at interior p points. periodicity handled
   %
   % @brief Method for vortex visualization
   %
   %
   % @details
   %
   % Method to compute the second invariant of the velocity
   % gradient tensor (Q), used to visualize vortices in 
   % complex and turbulent flows
   % 
   % Only valid for uniform cartesian meshes.
   % 
   % @code
   % fr = Frame3D('test1','path','/home/oflores/data');
   % fr.secondinvariant
   % 
   % isosurface(fr.yp,fr.xp,fr.zp,fr.Q,qref)
   %
   % @endcode
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   overrides = true;
   misc.assigndefaults(varargin{:});
   self.readfields('items',{'params','mesheu','buffers','bc'},...
                            'overrides',overrides);
   
   Q = zeros(self.sizebyvar.p.count(1)-2,...
             self.sizebyvar.p.count(2)-2,...
             self.sizebyvar.p.count(3)-2);

   perx = 0; 
   pery = 0; 
   perz = 0; 
   if self.bcux11 == 5
      perx = 1;
   elseif self.bcux21 == 5
      pery = 1;
   elseif self.bcux31 == 5
      perz = 1;
   end

   if ~self.isref
      dx = self.xux(2)-self.xux(1);
      dy = self.yux(2)-self.yux(1);
      dz = self.zux(2)-self.zux(1);
   else
      dx_init = diff(self.xux); dx = repmat(dx_init,[1,self.ny,self.nz]);
      dy_init = diff(self.yuy); dy = repmat(dy_init,[1,self.nx,self.nz]);
      dy = permute(dy,[2 1 3]);
      dz_init = diff(self.zuz); dz = repmat(dz_init,[1,self.ny,self.nx]);
      dz = permute(dz,[3 2 1]);
   end

   Q = -0.5*((diff(self.ux(1:end,2:end-1,2:end-1),1,1)./dx).^2 +  ... 
             (diff(self.uy(2:end-1,1:end,2:end-1),1,2)./dy).^2 +  ... 
             (diff(self.uz(2:end-1,2:end-1,1:end),1,3)./dz).^2 );
   
   
   [xu,yu,zu] = ndgrid(self.xux(1:end-perx),self.yux,self.zux);
   [xv,yv,zv] = ndgrid(self.xuy,self.yuy(1:end-pery),self.zuy);
   [xw,yw,zw] = ndgrid(self.xuz,self.yuz,self.zuz(1:end-perz));
   
   % dvdx * dudy: 
   % v to xux,yux,zux
   vv = interpn(self.xuy,self.yuy(1:end-pery),self.zuy,...
                self.uy(:,1:end-pery,:),xu,yu,zu);  
   % u to xuy,yuy,zuy
   uu = interpn(self.xux(1:end-perx),self.yux,self.zux,...
                self.ux(1:end-perx,:,:),xv,yv,zv);  
   
   Q = Q - ...
       diff(vv(:,2:end-1,2:end-1),1,1)./dx.*diff(uu(2:end-1,:,2:end-1),1,2)./dy;
   
   
   % dwdx * dudz: 
   % w to xux,yux,zux
   ww = interpn(self.xuz,self.yuz,self.zuz(1:end-perz),...
        self.uz(:,:,1:end-perz),xu,yu,zu);
   % u to xuz,yuz,zuz
   uu = interpn(self.xux(1:end-perx),self.yux,self.zux,...
                self.ux(1:end-perx,:,:),xw,yw,zw);
   
   Q = Q - ...
     diff(ww(:,2:end-1,2:end-1),1,1)./dx.*diff(uu(2:end-1,2:end-1,:),1,3)./dz;
   
   % dwdy * dvdz: 
   % w to xuy,yuy,zuy
   ww = interpn(self.xuz,self.yuz,self.zuz(1:end-perz),...
                self.uz(:,:,1:end-perz),xv,yv,zv);
   % v to xuz,yuz,zuz
   vv = interpn(self.xuy,self.yuy(1:end-pery),self.zuy,...
                self.uy(:,1:end-pery,:),xw,yw,zw);  
   
   Q = Q - ...
     diff(ww(2:end-1,:,2:end-1),1,2)./dy.*diff(vv(2:end-1,2:end-1,:),1,3)./dz;


   if isempty(find(strcmp(properties(self),'Q')))
      p = addprop(self,'Q');
   end
   self.Q = Q;
end 

function getlambda(self,varargin)
   %   @author M.Moriche
   %   @date 13-01-2014 by M.Moriche\n
   %         created
   %   @brief  lambda2 criterion calculation
   %
   %   @details
   %
   %   o  -  o  -  o  -  o  -  o  -  o      o/+  pressure points
   %                                        - ux points
   %   |     |     |     |     |     |      | uy points
   %                                  
   %   o  -  +  -  +  -  +  -  +  -  o      +  grad(u) 
   %                                       
   %   |     |     |     |     |     |
   %                                  
   %   o  -  +  -  +  -  +  -  +  -  o 
   %                                  
   %   |     |     |     |     |     |
   %                                  
   %   o  -  o  -  o  -  o  -  o  -  o 
   %
   %         +-                  -+
   %         |  dudx, dudy, dudz  |
   % gradu = |                    |
   %         |  dvdx, dvdy, dvdz  |
   %         |                    |
   %         |  dwdx, dwdy, dwdz  |
   %         +-                  -+
   %
   % x = s^2+o^2 = [a b; c d];
   % det(x-lambda*i) = 0
   %
   % lambda = ((a+d)-sqrt( (a+d).^2 + 4*(b.*c)))/2
   %
   % optional arguments: 
   %   - overrides: boolean to overrid existing fields.
   %
   %  usage example:
   %
   %  @code
   %  fr = frame2d(basenm, 'path', datapath);
   %  fr.getlambda();
   %  uu = mean(fr.lambda(:)) - std(fr.lambda(:));
   %  contour(fr.xl2, fr.yl2, fr.l2', [uu uu],'color','k')
   %  @endcode
   %
   
   
   overrides=true;
   misc.assigndefaults(varargin{:});
   
   self.readfields('items',{'mesheu','ux','uy','uz'},'overrides',overrides);
   
   nx = length(self.xp);
   ny = length(self.yp);
   nz = length(self.zp);
   ipos = 2:nx-1;
   jpos = 2:ny-1;
   kpos = 2:nz-1;
   nx = length(ipos);
   ny = length(jpos);
   nz = length(kpos);
   
   l2  = zeros(nx,ny,nz);
   xl2 = self.xp(ipos);
   yl2 = self.yp(jpos);
   zl2 = self.zp(kpos);
   
   dudx = diff(self.ux, 1, 1)/self.dx;
   % average in the direction of stagger
   uxprom = 0.5*self.ux(1:end-1,:,:) + 0.5*self.ux(2:end,:,:);
   % average in the direction of derivative 
   uxprom = 0.5*uxprom(:,1:end-1,:) + 0.5*uxprom(:,2:end,:);
   dudy = diff(uxprom, 1, 2)/self.dy;
   % average in the direction of stagger
   uxprom = 0.5*self.ux(1:end-1,:,:) + 0.5*self.ux(2:end,:,:);
   % average in the direction of derivative 
   uxprom = 0.5*uxprom(:,:,1:end-1) + 0.5*uxprom(:,:,2:end);
   dudz = diff(uxprom, 1, 3)/self.dz;
   clear uxprom
   
   % average in the direction of stagger
   uyprom = 0.5*self.uy(:,1:end-1,:) + 0.5*self.uy(:,2:end,:);
   % average in the direction of derivative 
   uyprom = 0.5*uyprom(1:end-1,:,:) + 0.5*uyprom(2:end,:,:);
   dvdx = diff(uyprom, 1, 1)/self.dx;
   dvdy = diff(self.uy, 1, 2)/self.dy;
   % average in the direction of stagger
   uyprom = 0.5*self.uy(:,1:end-1,:) + 0.5*self.uy(:,2:end,:);
   % average in the direction of derivative 
   uyprom = 0.5*uyprom(:,:,1:end-1) + 0.5*uyprom(:,:,2:end);
   dvdz = diff(self.uy, 1, 3)/self.dz;
   
   % average in the direction of stagger
   uzprom = 0.5*self.uz(:,:,1:end-1) + 0.5*self.uz(:,:,2:end,:);
   % average in the direction of derivative 
   uzprom = 0.5*uzprom(1:end-1,:,:) + 0.5*uzprom(2:end,:,:);
   dwdx = diff(uzprom, 1, 1)/self.dx;
   % average in the direction of stagger
   uzprom = 0.5*self.uz(:,:,1:end-1) + 0.5*self.uz(:,:,2:end,:);
   % average in the direction of derivative 
   uzprom = 0.5*uzprom(:,1:end-1,:) + 0.5*uzprom(:,2:end,:);
   dwdy = diff(uzprom, 1, 2)/self.dy;
   dwdz = diff(self.uz, 1, 3)/self.dz;
   
   
   for k = 1:nz
      for j = 1:ny
         for i = 1:nx
            a = [dudx(i,j,k),dudy(i,j,k),dudz(i,j,k);...
                 dvdx(i,j,k),dvdy(i,j,k),dvdz(i,j,k);...
                 dwdx(i,j,k),dwdy(i,j,k),dwdz(i,j,k)];
            s = 0.5*(a+transpose(a));
            o = 0.5*(a-transpose(a));
            x = s^2 + o^2;
            l = sort(eig(x));
            l2(i,j,k) = real(l(2));
         end
      end
   end
    
   sets = {'l2','xl2','yl2','zl2'};
   for i = 1:length(sets)
      if isempty(find(strcmp(properties(self),sets{i})))
         p = addprop(self,sets{i});
         self.(sets{i}) = misc.getval(sets{i});
      elseif overrides
         self.(sets{i}) = misc.getval(sets{i});
      end
   end
   
   
   
   
   
   return
end


end  % methods
end  % class

