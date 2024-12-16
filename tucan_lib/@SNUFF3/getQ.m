function getQ(self, varargin) 
   % @author O. FLores
   % @date 03-07-2014 by O. Flores \n 
   %       Created
   % @date 13-12-2014 by M. Moriche \n 
   %       Modified. Q at interior p points. periodicity handled
   % @date Modified 06-05-2016 by A.Gonzalo \n
   %       Adapt function to Non-Uniform cartesian meshes
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
   % @code
   % fr = Frame3D('test1','path','/home/oflores/data');
   % fr.secondinvariant
   % 
   % isosurface(fr.yp,fr.xp,fr.zp,fr.Q,qref)
   %
   % @endcode
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   Q = zeros(self.sizebyvar.p.count(1)-2,...
             self.sizebyvar.p.count(2)-2,...
             self.sizebyvar.p.count(3)-2);

   dcx = diff(self.xux);
   dcy = diff(self.yuy);
   dcz = diff(self.zuz);

   [dx dy dz] = ndgrid(dcx,dcy,dcz);

   perx = 0; 
   pery = 0; 
   perz = 0; 

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

   self.Q = Q;
end 

