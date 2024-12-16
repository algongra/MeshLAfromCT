classdef Run < dynamicprops

% \date 05-04-2015 by M.Moriche \n
%       Modified Frame. Added extension as optional input

properties
   path
   basenm
   frames
   raws
   stats
   ndim
end

properties(Access='protected')
coords = {'x','y','z'}
end

methods
function self = Run(basenm, ndim, varargin)
% \date 05-04-2015 by M.Moriche \n
%       Modified Frame. Added extension as optional input
   
   path = '.';
   frame_ext = 'h5frame';
   misc.assigndefaults(varargin{:});
   self.frames = NumberedFile(path,basenm,frame_ext, 2);
   self.raws   = NumberedFile(path,basenm,'h5raw',   2);
   self.stats  = NumberedFile(path,basenm,'h5stats',   2);

   self.path = path;
   self.basenm = basenm;
   self.ndim = ndim;
 
return
end

function ifrlist = getifr(self,item)
   % @author M.Moriche
   % @date Modified 25-05-2013 by M.Moriche
   % @brief get available frames saved
   %
   % @details 
   %
   % Reads the directory self.path and find the frames available.
   % making use of the generic function getifr.
   %
   % Uses regular expresion regfrmt and regfrmt_ifrpos to identify the
   % file.
   %
   % Example: 
   %
   % @code
   % fr = Frame(basename);
   % ifrlist = fr.getifr();
   % @endcode
   %
   dirnm = fullfile(self.path,self.basenm);
   ifrlist = misc.getifr(dirnm, self.(item).regex_ifr, self.(item).ifrpos);
end

function ifr = get(self, item, ifr)

   if isstr(ifr)
      ifrlist = self.getifr(item);
      if strcmp(ifr,'last')
         ifr = ifrlist(end);
      elseif strcmp(ifr,'first')
         ifr = ifrlist(1);
      end
   end
    
return
end

function tmean =  unitarytime2(self,varargin)
   % method to get 
   %   time/point/step*nproc
   %   @date 23-11-2014 by M.Moriche
   %         Bug in indices
   
   ifrlist = self.getifr('frames');
   misc.assigndefaults(varargin{:});
   
   items = {'istep'};
   for i = 1:self.ndim
      myvar1 = ['n' self.coords{i}];
      myvar2 = ['nd' self.coords{i}];
      items = [items myvar1 myvar2];
   end
   
   timesum  = 0.;
   time2sum = 0.;
   tmean = zeros(1, length(ifrlist) - 1);
   ii = 1;
   ifrlist
   for i1 = 1:length(ifrlist)-1
      ifr = ifrlist(i1);
      fr1 = Frame(self.basenm, self.ndim, 'path',self.path, 'ifr', ifr);
      fr1.readfields('items',items);
      step1 = double(fr1.istep);
   
      
      ifr = ifrlist(i1+1);
      fr2 = Frame(self.basenm, self.ndim, 'path',self.path, 'ifr', ifr);
      fr2.readfields('items',items);
      step2 = double(fr2.istep);
      nsteps = step2 - step1;
   
      npoints = 1;
      nprocs  = 1;
      for i2 = 1:self.ndim
         npoints = npoints*double(getfield(fr2,['n' self.coords{i2}]));
         nprocs  = nprocs*double(getfield(fr2,['nd' self.coords{i2}]));
      end
   
      ifr = ifrlist(i1);
      fnm = sprintf(self.frames.frmt, self.basenm ,ifr, 'h5frame');
      fnm = fullfile(self.path, self.basenm, fnm);
      uu = dir(fnm);
      t0 = datevec(uu.date);
   
      ifr = ifrlist(i1+1);
      fnm = sprintf(self.frames.frmt, self.basenm ,ifr, 'h5frame');
      fnm = fullfile(self.path, self.basenm, fnm);
      uu = dir(fnm);
      t1 = datevec(uu.date);
      dt = etime(t1,t0);
      dt = dt/(nsteps*npoints)*nprocs;
   
      tmean(ii) = dt;
      ii = ii + 1;
   end
   
   return
end

function bytes = getdiskspace(self,varargin)
   % @author M.Moriche
   % @date 21-10-2014
   % @brief Method to get the disk space in bytes used by the simulation
   %
   % @details
   %
   % @code
   % fr = Frame2D(basenm, 'path', datapath);
   % bb = fr.getdiskspace();   % bytes
   % Gb = bb/(2^30);           % Gyga bytes
   % @endcode
  
   item = 'all';
   misc.assigndefaults(varargin{:}); 

   if strcmp(item,'all')
      itemlist = {'frames','raws','stats'};
   else
      itemlist = {item};
   end
   
   mypatt = '(.*\n)*(\d*)(.*)total';

   bytes = 0;
   for i1 = 1:length(itemlist)
      item = itemlist{i1};
      [st,rs] = system(sprintf('du -sb --total %s/%s/*%s', self.path, self.basenm,self.(item).ext));
      if st == 0
         a = regexp(rs, mypatt, 'tokens');
         bb = str2num(a{1}{2});
         bytes = bytes + bb;
      else
         display('error... ')
         display(rs)
      end
   end
end

end %methods

end %classdef
