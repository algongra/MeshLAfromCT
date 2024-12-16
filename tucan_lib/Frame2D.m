classdef Frame2D < dynamicprops & Frame
%
% @author M.Moriche
%
% @brief Class to ease the manage of 2D output frames from TUCAN
% @date 08-09-2014 by M.Moriche\n
%       Added getlambda
%
% @details
%
% METHODS:
% - plotdomain()
% - streamplot()
% - vorticity()
% - stress()
% - collocate()
% - contourplot()
% - getlambda()
%
% EXAMPLES:
% =========
%
% Create a contour plot for the velocity:
% ---------------------------------------
%
% This example crates a contour plot of the fields previously read
%
% @code
% basenm = 'MASF_1004';
% datapath = '/data2/mmoriche/IBcode/MAS';
% fr = Frame2D(basenm, 'path', datapath);
% fr.readfields('items',{'time','xux','yux','ux'});
% pcolor(fr.xux, fr.yux, fr.ux');
% @endcode
%
% Create a contour plot for the velocity of part of the domain:
% -------------------------------------------------------------
%
% This example crates a contour plot with a own method and only
%  for part of the domain
%
% @code
% basenm = 'MASF_1004';
% datapath = '/data2/mmoriche/IBcode/MAS';
% slab = [-1, 5; -3, 3]; % [x0, xf; y0, yf];
% fr = Frame2D(basenm, 'path', datapath, 'slab', slab);
% fr.contourplot('wz','vmin',-5, 'vmax', 5);
% @endcode
properties
end
properties(Access='protected')
end
methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function self = Frame2D(basenm,varargin)

   ndim = 2;
   self@Frame(basenm,ndim,varargin{:});

end

%________________________________________________________________________________
%________PLOT DOMAIN_____________________________________________________________
%________________________________________________________________________________

function fig = plotdomain(self,varargin)
   %    \author M.Moriche
   %    \date Modified 25-05-2013 by M.Moriche
   %    \date 18-06-2014 by M.Moriche \n
   %          Modified
   %    \brief Plots an outline of the simulation domain
   %
   %  \details
   %  
   % No mandatory arguments
   %
   % Optional arguments:
   %  - fig: figure handle of the figure to be plotted  
   %  - widthratio: Ratio of the figure width wrt to the
   %     default value (docwidth=16.56/2.54)\
   %    W = docwidth*widthratio;
   %
   %  Usage example: 
   %
   %  \code
   %  fr = Frame2D(basename);
   %  fig = fr.plotdomain();
   %  \endcode
   %
   fig = [];
   widthratio=1/3;
   visible = 'on';
   misc.assigndefaults(varargin{:});
   if isempty(fig)
      fig = figure('visible',visible);
   else
      set(0,'CurrentFigure',fig);
   end
     
   needs = {'x0','y0','xf','yf','xib','yib'};
   self.readfields('items',needs,'overrides',true);

   hold on
   plot(self.xib,self.yib);
   axis equal
   axis off

   x = [self.x0,self.xf,self.xf,self.x0,self.x0];
   y = [self.y0,self.y0,self.yf,self.yf,self.y0];
   plot(x,y);
   xlim([self.x0 self.xf]);
   ylim([self.y0 self.yf]);

   lax = min(self.yf-self.y0,self.xf-self.x0)/10.;
   myplot.arrow([0,0],[lax,0])
   myplot.arrow([0,0],[0,lax])
   text(lax,0,'X')
   text(0,lax,'Y')
   rhb = (self.yf-self.y0)/(self.xf-self.x0);
   docwidth=16.56/2.54;
   W = docwidth*widthratio;
   % adjust figure size
   set(fig,'Units','inches');
   set(fig,'Position',[0 0 W W*rhb]);
   ax = gca;
   set(ax,'Units','Normalized');
   set(ax,'Position',[0 0 1 1]);
return
end

function fig = streamplot(self,varargin)
   %    \author M.Moriche
   %    \date Modified 25-05-2013 by M.Moriche
   %    \brief Plot streamlines
   %
   %   \details
   %
   % No mandatory arguments
   %
   % Optional arguments: 
   %  - fig: figure handle of the figure to use
   %  - showib: boolean to show or not the immersed
   %            boundary (DEFAULT: true)
   %  - mask: boolean to plot or not streamlines
   %              inside the ib (DEFAULT: false)
   %  - density: number of lines per side of box.
   %  - color: rgb array for the streamlines color
   %
   %  Usage example:
   %
   %  \code
   %  fr = Frame(basename)
   %  fig = fr.streamplot();
   %  fig = fr.streamplot("showib",false);
   %  fig = fr.streamplot('mask',true);
   %  fig = fr.streamplot('mask',true,'density',10);
   %  fig = fr.streamplot('color',[1 0 0]);
   % \endcode
   %

   fig = [];
   showib=true;
   mask=true;
   density=30;
   color=[0.5 0.5 0.5]; % rgb 
   overrides=true;
   visible='on';
   misc.assigndefaults(varargin{:});
   if showib | mask
      needs = {'xib','yib'};
      self.readfields('items',needs,'overrides',overrides);
   end

   % xpol, ypol define the polygon within streamlines
   %  are not to be started.
   if mask
      xpol = self.xib;
      ypol = self.yib;
   else
      %xpol = [0,0.00001,0];
      %ypol = [0,0,0.00001];
      %xpol = [NaN,NaN,NaN];
      %ypol = [NaN,NaN,NaN];
      xpol = [0,1,0];
      ypol = [0,0,1];
   end

   % needs to be collocated
   if isempty(find(strcmp(properties(self),'xcoll')))
      self.collocate('overrides',true);
   end
         
   if isempty(fig)
      fig = figure('visible',visible);
   else
      set(0,'CurrentFigure',fig);
   end

   [Y,X]=meshgrid(self.ycoll,self.xcoll);
   mingap = 0.02*min([self.xcoll(end)-self.xcoll(1),self.ycoll(end)-self.ycoll(1)]);

   xiarr = linspace(self.xcoll(1),self.xcoll(end),density);
   yiarr = linspace(self.ycoll(1),self.ycoll(end),density);

   xacc = [999999];
   yacc = [999999];
   for i=1:length(xiarr) 
      xi = xiarr(i);
      for j=1:length(yiarr) 
         yi = yiarr(j);
         in = inpolygon(xi,yi,xpol,ypol);
         if ~in
            dd = sqrt( (xi-xacc).^2 + (yi-yacc).^2);
            if min(dd) > mingap
               s1=streamline(X',Y', self.uxcoll', self.uycoll',xi,yi);
               s2=streamline(X',Y',-self.uxcoll',-self.uycoll',xi,yi);
               set(s1,'Color',color);
               set(s2,'Color',color);

               xacc = [xacc get(s1,'XData')];
               xacc = [xacc get(s2,'XData')];
               yacc = [yacc get(s1,'YData')];
               yacc = [yacc get(s2,'YData')];
            end
         end
      end
   end

   %axis equal
   %xlabel('x')
   %ylabel('y')
end

function vorticity(self)
   %%   @author M.Moriche
   %    @dated 22-07-2017 by M.Moriche\n
   %           Created
   %    @brief Vorticity calculation
   %
   %    @details
   %
   %    Uses inheritance from super class Frame
  
   vorticity@Frame(self, 3);

return
end
function stress(self,varargin)
   %%   \author M.Moriche
   %    \date Modified 25-05-2013 by M.Moriche
   %    \brief Calculate viscous stress tensor (Eulerian frame)
   %
   %  \details
   %
   %  \code
   %  o  -  o  -  o  -  o  -  o  -  o      o/+  pressure points
   %                                       - ux points
   %  |  *  |  *  |  *  |  *  |  *  |      | uy points
   %                                 
   %  o  -  +  -  +  -  +  -  +  -  o      + tauxx,tauyy (interior p)
   %                                       * tauxy (x of ux, y of uy)
   %  |  *  |  *  |  *  |  *  |  *  |
   %                                 
   %  o  -  +  -  +  -  +  -  +  -  o 
   %                                 
   %  |  *  |  *  |  *  |  *  |  *  |
   %                                 
   %  o  -  o  -  o  -  o  -  o  -  o 
   %
   %  \endcode
   %
   %  NOTE: tauxx,tauyy are calculated in p points while
   %     tauxy is in vorticity points (x of ux, y of uy)
   %  
   %  Usage example:
   %
   %  \code
   %  fr = Frame(basename);
   %  fr.stress()
   %  \endcode
   %
   overrides = true;
   misc.assigndefaults(varargin{:});

   % read fields that it needs from current ifr
   needs = {'xux','yuy','ux','uy','dx','dy','Re'};
   self.readfields('items',needs,'overrides',overrides);
   % add dynamic props
   sets = {'tauxx','tauxy','tauyy'};
   for i = 1:length(sets)
      if isempty(find(strcmp(properties(self),sets{i})))
         p = addprop(self,sets{i});
      end
   end

   tauxx= 2.*(self.ux(2:end,2:end-1) - self.ux(1:end-1,2:end-1))/self.dx;
   tauyy= 2.*(self.uy(2:end-1,2:end) - self.uy(2:end-1,1:end-1))/self.dy;
   tauxy= (self.uy(2:end,:) - self.uy(1:end-1,:))/self.dx...
        + (self.ux(:,2:end) - self.ux(:,1:end-1))/self.dy;
      
   self.tauxx = tauxx/self.Re;
   self.tauyy = tauyy/self.Re;
   self.tauxy = tauxy/self.Re;
end

function collocate(self,varargin)
   %    \author M.Moriche
   %    \date Modified 25-05-2013 by M.Moriche
   %    \brief   Interpolates velocities to have collocated mesh data
   %  
   %  \details
   %
   %  \code
   %   o  -  o  -  o  -  o  -  o  -  o      o/+  pressure points
   %                                        - ux points
   %   |     |     |     |     |     |      | uy points
   %                                  
   %   o  -  +  -  +  -  +  -  +  -  o      +  (p,ux,uy) collocated
   %                                       
   %   |     |     |     |     |     |
   %                                  
   %   o  -  +  -  +  -  +  -  +  -  o 
   %                                  
   %   |     |     |     |     |     |
   %                                  
   %   o  -  o  -  o  -  o  -  o  -  o 
   %  
   %  \endcode
   %
   %   Usage example: 
   % 
   %  \code
   %   fr = Frame(basename);
   %   fr.collocate();
   %   quiver(fr.xcoll,fr.ycoll,fr.uxcoll,fr.uycoll)
   %  \endcode
   %

   overrides=true;
   misc.assigndefaults(varargin{:});
   % read fields that it needs from current ifr
   needs = {'xp','yp','ux','uy','p'};
   self.readfields('items',needs,'overrides',overrides);
   % add dynamic props
   sets = {'xcoll','ycoll','uxcoll','uycoll','pcoll'};

   uxcoll = 0.5*(self.ux(1:end-1,2:end-1) + self.ux(2:end,2:end-1));
   uycoll = 0.5*(self.uy(2:end-1,1:end-1) + self.uy(2:end-1,2:end));
   pcoll = self.p(2:end-1,2:end-1);
   xcoll = self.xp(2:end-1);
   ycoll = self.yp(2:end-1);

   for i = 1:length(sets)
      if isempty(find(strcmp(properties(self),sets{i})))
         p = addprop(self,sets{i});
         self.(sets{i}) = misc.getval(sets{i});
      elseif overrides
         self.(sets{i}) = misc.getval(sets{i});
      end
   end
   self.uxcoll = uxcoll;
   self.uycoll = uycoll;
   self.pcoll = pcoll;
   self.xcoll = xcoll;
   self.ycoll = ycoll;
end

function fig = contourplot(self,var,varargin)
   %   @author M.Moriche
   %   @date Modified 25-05-2013 by M.Moriche
   %   @brief  Contour plot of specific variables
   %
   %   @details
   %
   % Mandatory arguments
   %  - variable: p,ux,uy,w
   %
   % Optional arguments: 
   %  - fig: figure handle of the figure to use
   %  - showib: boolean to show or not the immersed
   %            boundary (DEFAULT: true)
   %  - mask: boolean to plot or not
   %              inside the ib (DEFAULT: true)
   %  - vmin: minimum value for the caxis (DEFAULT: Min of field)
   %  - vmax: maximum value for the caxis (DEFAULT: Max of field)
   %  - fast: boolean to do the contour with imagesc or pcolor
   %  - colorbarwidth: float setting the colorbar widht to a ratio
   %     with respect to the figure axes width.
   %     DEFAULT: 0.05 (set to 0 if nothing is to be done)
   %
   %  Usage example:
   %
   %  @code
   %  fr = Frame2D(basename)
   %  fig = fr.contourplot('p');
   %  fig = fr.contourplot('p''showib',false);
   %  fig = fr.contourplot('ux''mask',true);
   %  @endcode
   %
   fig = [];
   vmin = [];
   vmax = [];
   fast = true;
   showib = true;
   mask = true;
   visible = 'on';
   colorbarwidth=0.05;
   misc.assigndefaults(varargin{:});
   
   if isempty(fig)
      fig = figure('visible',visible);
   else
      set(0,'CurrentFigure',fig);
   end
   hold on;
   
   if showib | mask
      needs = {'xib','yib','body_ib','body_ie'};
      self.readfields('items',needs,'overrides',false);
   end
   
   % xpol, ypol define the polygon within streamlines
   %  are not to be started.
   if mask
      xpol = self.xib;
      ypol = self.yib;
   else
      xpol = [0,0.000001,0];
      ypol = [0,0,0.0000001];
   end
   
   if strcmp(var,'wz')
      self.vorticity();
      self.readfields('items',{'xux','yuy'},'overrides',false);
      x = self.xux;
      y = self.yuy;
   else
      needs = {var,['x' var],['y' var],'dx','dy'};
      self.readfields('items',needs,'overrides',false);
      x = self.(['x' var]);
      y = self.(['y' var]);
   end
   self.readfields('items',{'x0','xf','y0','yf'},'overrides',false);
   if fast
      x2 = x;
      y2 = y;
   else
      x2 = [x-self.dx/2;x(end)+self.dx/2];
      y2 = [y-self.dy/2;y(end)+self.dy/2];
   end
   
   buff = self.(var);
   [Y X] = meshgrid(y,x);
   in = inpolygon(X,Y,xpol,ypol);
   buff(in) = NaN;
   
   if isempty(vmin)
      vmin = min(buff(:));
   end
   if isempty(vmax)
      vmax = max(buff(:));
   end
   colormap(myplot.redblue);
   if fast
      imagesc(x2,y2,buff');
   else
      buff = [buff ones(size(buff,1),1)];
      buff = [buff;ones(1,size(buff,2))];
      h = pcolor(x2,y2,buff');
      set(h,'edgecolor','none');
   end
   
   caxis([vmin vmax]);
   set(gca,'YDir','normal');
   
   if showib
      try
         for i = 1:length(self.body_ib)
         
            ib = self.body_ib(i);
            ie = self.body_ie(i);
            fill(self.xib(ib:ie),self.yib(ib:ie),[0.4 0.4 0.4]);
         end
      catch
            fill(self.xib,self.yib,[0.4 0.4 0.4]);
      end
   end
   
   %xlim([self.x0 self.xf]);
   %ylim([self.y0 self.yf]);
   
   set(gca,'XTick',[self.x0 self.xf]);
   set(gca,'YTick',[self.y0 self.yf]);
   axis equal;
   axis tight
   c = colorbar;
   if colorbarwidth > 0
      xxx1  = get(gca,'position');
      xxx  =  get( c ,'position');
      xxx(3) = xxx1(3)*colorbarwidth;
      set(c,  'position',xxx);
      set(gca,'position',xxx1);
   end
   
   xlabel('x');
   ylabel('y');
   return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function getlambda(self,varargin)
   %   @author m.moriche
   %   @date 08-09-2014 by m.moriche\n
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
   %         +-            -+
   %         |  dudx, dudy  |
   % gradu = |              |
   %         |  dvdx, dvdy  |
   %         +-            -+
   %
   %       +-                          -+
   %       |               dudy   dvdx  |
   %       |      dudx,    ---- + ----  |
   %       |                2      2    |
   %  s =  |                            |
   %       |  dudy   dvdx               |
   %       |  ---- + ----,     dvdy     |
   %       |   2      2                 |
   %       +-                          -+
   %      
   %
   %
   %
   %       +-                          -+
   %       |               dudy   dvdx  |
   %       |       0,      ---- - ----  |
   %       |                2      2    |
   %  o =  |                            |
   %       |  dvdx   dudy               |
   %       |  ---- - ----,      0       |
   %       |   2      2                 |
   %       +-                          -+
   %      
   %               +-                                            -+
   %               |       2      2      2                        |
   %               |  - o12  + s11  + s12 ,   s11 s12 + s12 s22   |
   %   s^2 + o^2 = |                                              |
   %               |                             2      2      2  |
   %               |    s11 s12 + s12 s22,  - o12  + s12  + s22   |
   %               +-                                            -+
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
   
   self.readfields('items',{'dx'},'overrides',overrides);
   
   self.collocate('overrides', overrides);
   
   nx = length(self.xcoll);
   ny = length(self.ycoll);
   ipos = 2:nx-1;
   jpos = 2:ny-1;
   
   
   nx = length(ipos);
   ny = length(jpos);
   lambda  = zeros(nx,ny);
   lambda_eig  = zeros(nx,ny);
   xl2 = self.xcoll(ipos);
   yl2 = self.ycoll(jpos);
   
   %for i = ipos
      %for j = jpos
   dudx = (self.uxcoll(ipos+1,jpos)-self.uxcoll(ipos-1,jpos))/(2*self.dx);
   dudy = (self.uxcoll(ipos,jpos+1)-self.uxcoll(ipos,jpos-1))/(2*self.dx);
   dvdx = (self.uycoll(ipos+1,jpos)-self.uycoll(ipos-1,jpos))/(2*self.dx);
   dvdy = (self.uycoll(ipos,jpos+1)-self.uycoll(ipos,jpos-1))/(2*self.dx);
   %for i = 1:nx
   %   for j = 1:ny
   %      a = [dudx(i,j),dudy(i,j);dvdx(i,j),dvdy(i,j)];
   %      s = 0.5*(a+transpose(a));
   %      o = 0.5*(a-transpose(a));
   %      x = s^2 + o^2;
   %      l = sort(eig(x));
   %      lambda_eig(i,j) = l(1);
   %
   %      %%%
   %      s11 = dudx(i,j);
   %      s12 = 0.5*(dudy(i,j)+dvdx(i,j));
   %      s22 = dvdy(i,j);
   %      o12 = 0.5*(dudy(i,j)-dvdx(i,j));
   %      a = s11.^2 + s12.^2 - o12.^2;
   %      b = s11.*s12 + s12.*s22;
   %      c = b;
   %      d = s12.^2 + s22.^2 - o12.^2;
   %      x2 = [a b; c d];
   %      lambda(i,j) = real(((a+d)-sqrt((a+d).^2 - 4*(a.*d - b.*c))  )/2;
   %
   %      %if lambda_eig(i,j) ~= lambda(i,j)
   %      %   lambda(i,j)
   %      %   lambda_eig(i,j)
   %      %   pause
   %      %end
   %      if ~ isreal(lambda(i,j))
   %         lambda(i,j)
   %         lambda_eig(i,j)
   %         pause
   %      end
   %   end
   %   i
   %end
    
   s11 = dudx;
   s12 = 0.5*(dudy+dvdx);
   s22 = dvdy;
   o12 = 0.5*(dudy-dvdx);
   a = s11.^2 + s12.^2 - o12.^2;
   b = s11.*s12 + s12.*s22;
   c = b;
   d = s12.^2 + s22.^2 - o12.^2;
   l2 = real(((a+d)-sqrt((a+d).^2 - 4*(a.*d - b.*c))  )/2);
   
   sets = {'l2','xl2','yl2'};
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

