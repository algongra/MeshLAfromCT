classdef Frame <  NumberedFile & dynamicprops
%
% \author M.Moriche
% \date 25-05-2013 by M.Moriche \n
%       Modified
% \date 07-06-2013 by M.Moriche\n
%       Modified 
% \date 07-06-2013 by M.Moriche \n
%       Imported from trunk/ibmpi and set to 
%       work with ibsharp. Basic change, add readh5frame
% \date 17-12-2013 by M.Moriche\n
%       Documented
% \date 16-05-2014 by M.Moriche\n
%       Modifications including Classes inheritance 
% \date 10-11-2014 by M.Moriche\n
%       Removed unitarytime2, is in Run methods \n
%       Modified unitarytime
% \date 24-12-2014 by M.Moriche \n
%       Modified Frame. added slab 
% \date 25-12-2014 by M.Moriche \n
%       Modified Frame. slab working.
% \date 05-04-2015 by M.Moriche \n
%       Modified Frame. Added extension as optional input
% \date 04-11-2015 by A.Gonzalo \n
%       Modified Frame. Added isref property
%
% \brief Class for simulation frames handling
%
%
% \details
%
% Class to simplify the postprocessing of simulation frames.
% Scripts that use this class result much more clear
%
% METHODS:
%  - readfields 
%  - setfield
%  - time = steptime
%  - tuni = unitarytime
%  - vorticity 
%  - bytes = getdiskspace
%  
%  
% \code
% datapath = '/data2/mmoriche/IBcode/...'
% slab = [-1, 2; -3, 3];
% fr = Frame('validation_stcyl',2);
% fr = Frame('validation_stcyl',2,'path', datapath);
% fr = Frame('validation_stcyl',2,'path', datapath, 'slab', slab);
% \endcode
%
properties
   sizebyvar
   isref
end
properties(Access='private')
   ndim;
end
properties(Access='protected')
   coords = {'x','y','z'};
   itemsByGroup;
end
methods
   function self = Frame(basenm,ndim,varargin)
      % @author M.Moriche
      % @date 17-12-2013 by M.Moriche \n
      %       Documented
      % @date 16-05-2014 by M.Moriche \n
      %       Renewed version with Class inheritance
      % @date 24-12-2014 by M.Moriche \n
      %       Added slab
      % @date 04-11-2015 by A.Gonzalo \n
      %       Added isref as an optional argument
      %       Changed old params (tsave,traw1,...) for new params (nstep,nraw,...)
      %
      % @brief Function to handle output frames from TUCAN
      %
      % @details
      %
      % MANDATORY ARGUMENTS:
      % - basenm: basename of the simulation
      % - ndim: rank of the problem (2 or 3)
      %
      % OPTIONAL ARGUMENTS:
      % - path = '.': directory where the simulation folder is placed
      % - ifr = 'last': index frame of the simulation to be read. 
      %                 Possible values are: integers, 'first', 'last'.
      % - slab = NaN(ndim,2); Box to make an slab of the fields
      % - doslab = true; Flag to slab the domain
      % - isref = false; Flag to know if the case has a non-uniform mesh
      % - ext = 'h5frame': extension of the file name
      % - ifrpos = 2: Position of the index in the groups from the SuperClass object
      %               NumberedFile.regex (Used to read availabel files)
      %
      % examples:
      %
      %  @code
      %  fr = Frame('validation_stcyl',2)
      %  fr = Frame('validation_stcyl',3,'ifr',12)
      %  fr = Frame('validation_stcyl',3,'path','/another/dir')
      %  fr = Frame('validation_stcyl',3,'ifr',23,'path','/another/dir')
      %  fr = Frame('validation_stcyl',3,'ifr','first')
      %  fr = Frame('validation_stcyl',2,'slab',[-1,3;-4,4],'ifr','first')
      %  fr = Frame('validation_stcyl',3,'slab',[-1,3;-4,4;-1.5,1.5],'ifr','first')
      %  @endcode
      %
   
      % optional arguments
      path = '.';
      ifr = 'last';
      slab = NaN(ndim, 2);
      doslab = true;
      isref = false;
      % setup for File Object
      ext = 'h5frame';
      ifrpos = 2;
      misc.assigndefaults(varargin{:});


      % SuperClass inheritance
      self@NumberedFile(path, basenm, ext, ifrpos);
   
      % achtung
      % get ifr
      r = Run(basenm, 2, varargin{:},'frame_ext',ext);
      ifr = r.get('frames',ifr);
      self.ifr  = ifr;
   
      % generate items to read fields (dependent on ndim)
      % BUFFERS
      buffers = {'p'};
      for i = 1:ndim
         myvar = ['u' self.coords{i}];
         buffers = [buffers myvar]; 
      end
      % MESHEU
      mesheu = {}; 
      for i = 1:ndim
         myvar = [self.coords{i} 'p'];
         mesheu = [mesheu myvar]; 
         for j = 1:ndim
            myvar = [self.coords{i} 'u' self.coords{j}];
            mesheu = [mesheu myvar]; 
         end
         if isref
            myvar = ['d' self.coords{i} 'r'];
         else
            myvar = ['d' self.coords{i}];
         end
         mesheu = [mesheu myvar]; 
      end
      % MESHLAG
      meshlag = {}; 
      for i = 1:ndim
         myvar = [self.coords{i} 'ib'];
         myvar2 = ['f' self.coords{i}];
         meshlag = [meshlag myvar myvar2]; 
      end
   
      % LAMBDA
      lambda = {'l2','l2c'};
      for i = 1:ndim
         myvar = [self.coords{i} 'l2'];
         lambda = [lambda myvar]; 
      end
   
      % bc
      bc = {};
      myvar = 'p';
      for idir = 1:ndim
         for ipos = 1:2
         myvar2 = ['bc' myvar num2str(idir) num2str(ipos)];
         bc = [bc myvar2];
         end
      end
      for ivar = 1:ndim
         myvar = ['u' self.coords{ivar}];
         for idir = 1:ndim
            for ipos = 1:2
            myvar2 = ['bc' myvar num2str(idir) num2str(ipos)];
            bc = [bc myvar2];
            end
         end
      end
      % PARAMS
      params = {'Re','nstep','geomfilename','nsdv','sdv',...
                'resfilename','svnversion','machine',...
                'makefile','ndim','nsave','nraw',...
                'nsamp','nstats'};
      for i = 1:ndim
          myvar1 = ['n'  self.coords{i}];
          myvar2 = ['nd' self.coords{i}];
          myvar3 = [self.coords{i} '0'];
          myvar4 = [self.coords{i} 'f'];
          myvar5 = ['bodyf' self.coords{i}];
          params = [params myvar1 myvar2 myvar3 myvar4];
      end
   
      itemsByGroup.mesheu = mesheu;
      itemsByGroup.meshlag = meshlag;
      itemsByGroup.params = params;
      itemsByGroup.buffers = buffers;
      itemsByGroup.lambda  = lambda ;
      itemsByGroup.bc  = bc ;
   
      self.ndim = ndim;
      self.itemsByGroup = itemsByGroup;
      self.isref = isref;
   
      % define box size for trunked variables
   
      self.sizebyvar = struct;
    
      if doslab
      self.readfields('items',{'mesheu'});
    
      % ij(ivar, idim, ipos)
      % ipos=1 beginning
      % ipos=2 ending
      ij = zeros(ndim+1,ndim,2); 
      ivar = 1; % set for p
      
      for ipos = 1:2
         for idim = 1:ndim
            xxname = [self.coords{idim} 'p'];
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
            ia = ia(1);
            ij(1,idim,ipos) = ia;
         end
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
         varnm = ['u' self.coords{ivar}];
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
            varnm = [self.coords{idim} 'u' self.coords{ivar}];
            sizebyvar.(varnm).offset = [ij(ivar+1,idim,1)-1];
            sizebyvar.(varnm).stride =                 [1];
            sizebyvar.(varnm).count =  [ij(ivar+1,idim,2)-ij(ivar+1,idim,1)+1];
            sizebyvar.(varnm).block =                  [1];
         end
      end
      for idim = 1:ndim
         varnm = [self.coords{idim} 'p'];
         sizebyvar.(varnm).offset = [ij(1,idim,1)-1];
         sizebyvar.(varnm).stride =                 [1];
         sizebyvar.(varnm).count =  [ij(1,idim,2)-ij(1,idim,1)+1];
         sizebyvar.(varnm).block =                  [1];
      end
      
      self.sizebyvar = sizebyvar;
   
      self.readfields('items',{'mesheu'});
      end
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   function readfields(self,varargin)
      % @author M.Moriche
      % @date Modified 25-05-2013 by M.Moriche
      % @brief  Reads fields of the frame
      %    
      % @details
      %
      % Do not know about global or local indices, therefore
      % global indices are taken and names modified to remove "_glob"
      % Also trimmed to remove spaces
      % 
      % NO MANDATORY ARGUMENTS
      %
      % OPTIONAL ARGUMENTS:
      % - items={'mesheu','meshlag','buffers'}: list of variables to be read.
      % - overrides: boolean to override an existing field (DEFAULT: true)
      %
      % @code
      % fr.readfields();
      % fr.readfields('items',{'xib','yib'});
      % fr.readfields('items',{'ux','uy'},'overrides',false);
      % fr.readfields('items',{'params','mesheu'},'overrides',false);
      % @endcode
      %
   
      % default items to read
      items = {'mesheu','buffers','meshlag'}; 
      % flag to override existing items
      overrides = true;
      misc.assigndefaults(varargin{:});
     
      groups = fields(self.itemsByGroup);
      items2 = {};
      for i = 1:length(items)
         fldnm = items{i};
         if ~ isempty(intersect(groups,{fldnm}))
            items2 = [items2 self.itemsByGroup.(fldnm)];
         else
            items2 = [items2 fldnm];
         end
      end
      items = items2;
   
      fnm = sprintf(self.frmt, self.basenm,self.ifr,self.ext);
      fnm = fullfile(self.path, self.basenm, fnm);
      [fr,failitems] = io.readh5file(fnm,items,'sizebyvar',self.sizebyvar);
      if length(failitems) > 0
         display('    items not found:')
         failitems
      end
   
      fldlist = fields(fr);
      for i =1:length(fldlist)
         fld = fldlist{i};
         if ~isempty(getfield(fr,fld))
            self.setfield(fld,getfield(fr,fld),'overrides',overrides);
         end
      end
   
   end 
   
   function setfield(self,fieldname,field,varargin)
      % @author M.Moriche
      % @date Modified 25-05-2013 by M.Moriche
      % @date 25-05-2013 by M.Moriche \n
      %       Modified
      % @brief  Sets a field manually (not read from data saved)
      %   
      % @details
      %
      %
      % Mandaroty arguments:
      %  - fieldname: character string
      %  - field:  field to be set
      %
      % Optinal arguments: 
      %  - overrides: boolean to override or not an existing field
      %     DEFAULT: true
      %
      % Usage example:
      %   
      % @code
      % fr = Frame(basename);
      % ux = average(cases.... blba bla bla);
      % fr.setfield('ux',ux);
      % fr.setfield('ux',ux,'overrides',false);
      % @endcode
      % 
   
      overrides=true;
      misc.assigndefaults(varargin{:});
   
      fieldname = strtrim(fieldname);
      fieldname = strrep(fieldname,'_glob','');
      if isempty(find(strcmp(properties(self),fieldname)))
         p = addprop(self,fieldname);
         self.(fieldname) = field;
      elseif overrides
         self.(fieldname) = field;
      end
   end

   function removefield(self, fieldname)
      m = findprop(self, fieldname);
      delete(m);
   return
   end

   function tstep =  steptime(self, varargin)
      % method to get 
      %   time/point/step*nproc
      %   @date 11-11-2014 by M.Moriche \n
      %         Added optional argument inidices
      %
      %   @details
      %
      %  NO MANDATORY ARGUMENTS
      %
      %  OPTIONAL ARGUMENTS:
      %  - indices = [1]: index list to sum the time contribution.
      
      indices = [1]; % default for total time
      misc.assigndefaults(varargin{:});
      
      items = {'istep_local','secvalues','secrate'};
      for i = 1:self.ndim
         myvar1 = ['n' self.coords{i}];
         myvar2 = ['nd' self.coords{i}];
         items = [items myvar1 myvar2];
      end
      
      self.readfields('items',items);
      
      nsteps  = double(self.istep_local+1);
      ttaux = sum(max(self.secvalues(indices,:)'));
      t = double(ttaux)/double(self.secrate);
      tstep = t/nsteps;
      
      return
   end

   function tuni =  unitarytime(self, varargin)
      % method to get 
      %   time/point/step*nproc
      %   @date 11-11-2014 by M.Moriche \n
      %         Added optional argument inidices
      %
      %   @details
      %
      %  NO MANDATORY ARGUMENTS
      %
      %  OPTIONAL ARGUMENTS:
      %  - indices = [1]: index list to sum the time contribution.
      
      indices = [1]; % default for total time
      misc.assigndefaults(varargin{:});
      
      items = {'istep_local','secvalues','secrate'};
      for i = 1:self.ndim
         myvar1 = ['n' self.coords{i}];
         myvar2 = ['nd' self.coords{i}];
         items = [items myvar1 myvar2];
      end
      
      self.readfields('items',items);
      
      npoints = 1;
      nprocs  = 1;
      for i = 1:self.ndim
         npoints = npoints*double(getfield(self,['n' self.coords{i}]));
         nprocs  = nprocs*double(getfield(self,['nd' self.coords{i}]));
      end
      nsteps  = double(self.istep_local+1);
      ttaux = sum(max(self.secvalues(indices,:)'));
      t = double(ttaux)/double(self.secrate);
      tuni = t/npoints/nsteps*nprocs;
      
      return
   end
   
   
   function vorticity(self,dim, varargin)
      % @author M.Moriche
      % @date Modified 25-05-2013 by M.Moriche
      % @brief  Gets vorticity field
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
      % fr = Frame(basename);
      % fr.vorticity();
      % imagesc(fr.w')
      % @endcode
      %
      overrides = true;
      misc.assigndefaults(varargin{:});
      dimlist = 1:self.ndim;
      dimstoread = setdiff(dimlist, dim);
      needs = {};
      for i = 1:length(dimstoread)
         nm = self.coords{dimstoread(i)};
         needs = [needs ['u' nm] ['d' nm]];
      end
      self.readfields('items',needs,'overrides',overrides);
   
      nm = self.coords{dim};
      dir1 = dimstoread(1);
      dir2 = dimstoread(2);
      nm1 = self.coords{dir1};
      nm2 = self.coords{dir2};
      fld1 = getfield(self, ['u' nm1]);
      fld2 = getfield(self, ['u' nm2]);
      h1   = getfield(self, ['d' nm1]);
      h2   = getfield(self, ['d' nm2]);
      %w = diff(self.uy,1,1)/self.dx - diff(self.ux,1,2)/self.dy;
      [i1,j1,k1] = size(diff(fld2,1,dir1));
      [i2,j2,k2] = size(diff(fld1,1,dir2));
      iend = min(i1,i2);
      jend = min(j1,j2);
      kend = min(k1,k2);
      d1 = diff(fld2,1,dir1)/h1;
      d2 = diff(fld1,1,dir2)/h2;
      w = d1(1:iend, 1:jend, 1:kend) - d2(1:iend, 1:jend, 1:kend);
      if isempty(find(strcmp(properties(self),['w' nm])))
         p = addprop(self,['w' nm]);
      end
      self.(['w' nm]) = w;
   
   end
   
   
   function bytes = getdiskspace(self)
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
   
      mypatt = '(\d*).*';
      [st,rs] = system(sprintf('du -sb %s', self.getfullfilename));
      if st == 0
         a = regexp(rs, mypatt, 'tokens');
         bytes = str2num(a{1}{1});
      else
         display('error... ')
         display(rs)
      end
   end

end  % methods
end  % class

