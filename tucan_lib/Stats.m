classdef Stats < dynamicprops
%   @author M.Moriche
%   @date Documented 09-05-2013 by M.Moriche
%   @date 22-08-2014 by M.Moriche \n
%         readfieldsh5
%   @brief Class to ease the handling of simualtion statistics
%
%  @details
%
% This class is used to read the statistics output files from 
%  TUCAN.
%
% Two extensions are handled:
% - h5stats
% - h5raw
%
% By default the class reads h5stats files, so set optional 
%  argument raw to true to read h5raw files.
%
% EXAMPLES:
% =========
%
% How to read statistics output from TUCAN:
% -----------------------------------------
%
% @code
% basenm = 'MASF_1003';
% datapath = '/data2/mmoriche/IBcode/MAS';
% st = Stats(basenm, 'path', datapath, 'raw', true);
% st.readfieldsh5('items', {'time', 'fysum'});
%
% fig = figure();
% plot(st.time, st.fysum);
%
% @endcode
%
% How to merge sequence of statistics files in one:
% -------------------------------------------------
%
% If the number of h5raw or h5stats files if big, 
% the task of reading fields might be time consuming, 
% the method merge is very useful in this situation.
%
% First, the wanted fields must be read.
%
% @code
% basenm = 'MASF_1003';
% datapath = '/data2/mmoriche/IBcode/MAS';
% st = Stats(basenm, 'path', datapath, 'raw', true);
% st.readfieldsh5('items', {'time', 'fysum'});
% @endcode
%
% Then the object is merged into a single file in the 
% simulation's directory. No index is put in the file, 
% just "basename.h5raw"
%
% @code
% st.merge();
% @endcode
%
% Then, in a different session, the data of this file
%  can be accessed, resulting in a much faster process
%  since the class reads a single file.
%
% @code
% basenm = 'MASF_1003';
% datapath = '/data2/mmoriche/IBcode/MAS';
% st = Stats(basenm, 'path', datapath, 'raw', true);
% st.readmerge('items', {'time', 'fysum'});
% @endcode
properties
   basenm
   ifrlist
   path
end
properties(Access='protected')
   filepattern
   fileformat
   mergeformat
   raw
end
methods 
function self = Stats(basenm, varargin)

   ifrlist = [];
   path = '.';
   raw=false;
   misc.assigndefaults(varargin{:});

   self.basenm = basenm;
   self.path = path;

   if raw
      patt='(\w*?)\.(\d+?)\.h5raw';
      mergefrmt='%s.h5raw';
      narray = misc.getformatlength(patt,'path',fullfile(path,basenm));
      frmt=['%s.%0' num2str(narray(2)) 'd.h5raw'];
   else
      patt='(\w*?)\.(\d+?)\.h5stats';
      mergefrmt='%s.h5stats';
      narray = misc.getformatlength(patt,'path',fullfile(path,basenm));
      frmt=['%s.%0' num2str(narray(2)) 'd.h5stats'];
   end

   self.filepattern = patt;
   self.fileformat  = frmt;
   self.mergeformat  = mergefrmt;

   if isempty(ifrlist)
      ifravailable = self.getifr();
      self.ifrlist = ifravailable;
   elseif strcmp(class(ifrlist),'double') && length(ifrlist) == 1
         ifravailable = self.getifr();
         if ifrlist > 0
            self.ifrlist = ifravailable(1:ifrlist);
         elseif ifrlist < 0
            dummy = length(ifravailable)+ifrlist;
            self.ifrlist = ifravailable(dummy:end);
         end
   elseif strcmp(class(ifrlist),'double') && length(ifrlist) == 2
         ifravailable = self.getifr()
         self.ifrlist = ifravailable(ifrlist(1):ifrlist(2));
   elseif strcmp(class(ifrlist),'double') && length(ifrlist)  > 2
         ifravailable = self.getifr();
         self.ifrlist = [];
         for ifr = ifravailable
            if ~isempty(find(ifrlist==ifr))
               self.ifrlist = [self.ifrlist ifr];
            end
         end
   end

   return 
end

function ifrlist = getifr(self)
   %get available frames saved
   basenm = self.basenm;
   patt = self.filepattern;
   ifrlist=[];
   flist = dir(fullfile(self.path,basenm));
   for i = 1:length(flist)
      f = flist(i).name;
      a = regexp(f,patt,'tokens');
      if length(a) > 0
         if length(a{1}) == 2
         bs = a{1}{1};
         pifr = str2num(a{1}{2});
         if strcmp(bs,basenm) && isempty(find(ifrlist-pifr==0))
            ifrlist = [ifrlist pifr];
         end
         end
      end
   end
   ifrlsit =  sort(ifrlist);
   return 
end


function readfieldsh5(self,varargin)
   %
   %
   % \date 22-08-2014 by M.Moriche \n
   %       Edited so NaN is set when the field is absent
   % \date 26-11-2014 by M.Moriche \n
   %       Added solver group
   % \date 26-11-2014 by M.Moriche \n
   %       Added forces group
   
   items0 = {'time','step'};
   items = {};
   misc.assigndefaults(varargin{:});
   
   if isempty(items)
      items = {'cl','cd','cs'};
   end
   items = [items0 items];
   
   itemsByGroup.solver = {...
     'itersux','itersuy','itersuz','itersphi', ...
     'finalresux','finalresuy','finalresuz','finalresphi'};
   
   itemsByGroup.forces = {...
     'fxsum','fysum','fzsum',...
     'mxsum','mysum','mzsum'};
   
   itemsByGroup.all = {'time', ...
   'step', ...
   'cfl', ...
   'cl', ...
   'cd', ...
   'cs', ...
   'of', ...
   'sc', ...
   'sumdiveru', ...
   'diffinout', ...
   'uxmax', ...
   'uymax', ...
   'uzmax', ...
   'itersux', ...
   'itersuy', ...
   'itersuz', ...
   'itersphi', ...
   'finalresux', ...
   'finalresuy', ...
   'finalresuz', ...
   'finalresphi', ...
   'sdvout', ...
   'prom_ux', ...
   'prom_uy', ...
   'prom_uz', ...
   'prom_uxux', ...
   'prom_uyuy', ...
   'prom_uzuz', ...
   'prom_uxuy', ...
   'prom_uxuz', ...
   'prom_uyuz'};
   
   groups = fields(itemsByGroup);
   
   items2 = {};
   for i = 1:length(items)
      fldnm = items{i};
      if ~ isempty(intersect(groups,{fldnm}))
         items2 = [items2 itemsByGroup.(fldnm)];
      else
         items2 = [items2 fldnm];
      end
   end
   items = items2;
   
   ifr = self.ifrlist(1);
   basenm = self.basenm;
   path = self.path;
   
   indx=zeros(1,length(items));
   notifr = [];
   for k = 1:length(self.ifrlist)
      ifr = self.ifrlist(k);
      fnm = sprintf(self.fileformat, basenm, ifr);
      fnm = fullfile(self.path, basenm, fnm);
      try
         [fr,failitems] = io.readh5file(fnm,items);
         indx(k) = misc.getlastzeroes(fr.time); 
      catch
         notifr = [notifr ifr];
      end
   end
   self.ifrlist = setdiff(self.ifrlist, notifr);
   
   for i = 1:length(items)
   
      ifr = self.ifrlist(1);
      fnm = sprintf(self.fileformat, basenm, ifr);
      fnm = fullfile(self.path, basenm, fnm);
      fr = io.readh5file(fnm,items);
   
      fld = items{i};
      %%%val = fr.(fld)';
      try 
         val = fr.(fld)';
      catch
         val = -NaN*fr.time';
      end
      %%%%
      if ndims(val) == 1
         val = val(1:indx(1));
      elseif ndims(val) == 2
         val = val(:,1:indx(1));
      end
      % add prop if it does not exist
      if isempty(find(strcmp(properties(self),fld)))
         p = addprop(self,fld);
      end
   
      for k = 2:length(self.ifrlist)
         ifr = self.ifrlist(k);
         fnm = sprintf(self.fileformat, basenm, ifr);
         fnm = fullfile(self.path, basenm, fnm);
         fr = io.readh5file(fnm,items);
         %%%%preval = fr.(fld)';
         try 
            preval = fr.(fld)';
         catch
            preval = NaN*fr.time';
         end
         %%%%
         if ndims(val) == 1
            val = [val preval(1:indx(k))]; 
         elseif ndims(val) == 2
            val = [val preval(:,1:indx(k))]; 
         end
      end 
      self.(fld) = val;
   end
   
   self
   for i = 1:length(items)
      fld = items{i};
      if length(find(isnan(self.(fld)))) == length(self.time)
         mp = findprop(self, fld);
         delete(mp);
      end
   end
   
   return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fnm = merge(self, varargin)
   % function to save stats class with read fields in to an unique file
   %      basename.h5raw
   % OPTIONAL ARGUMENT:
   %   fnm: file where the data is saved

   fnm = fullfile(self.path,self.basenm,sprintf(self.mergeformat,self.basenm));
   misc.assigndefaults(varargin{:});

   myfields = fields(self); 
   notsave = {'basenm','path'};
   
   tosave = setdiff(myfields,notsave);
   
   for i = 1:length(tosave) 
      varname = tosave{i};
      eval([varname '= [];']);
      misc.assignval(varname, self.(varname));
   end
   
   io.writeh5file(tosave,fnm);
   
   return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function readmerge(self,varargin)

   items = {};
   misc.assigndefaults(varargin{:});
   
   if isempty(items)
      items = {'all'};
   end
   itemsByGroup.all = {'ifrlist','time', ...
   'step', ...
   'cfl', ...
   'cl', ...
   'cd', ...
   'cs', ...
   'of', ...
   'sc', ...
   'sumdiveru', ...
   'diffinout', ...
   'uxmax', ...
   'uymax', ...
   'uzmax', ...
   'itersux', ...
   'itersuy', ...
   'itersuz', ...
   'itersphi', ...
   'finalresux', ...
   'finalresuy', ...
   'finalresuz', ...
   'finalresphi', ...
   'fxsum', ...
   'fysum', ...
   'mzsum', ...
   'sdvout', ...
   'prom_ux', ...
   'prom_uy', ...
   'prom_uz', ...
   'prom_uxux', ...
   'prom_uyuy', ...
   'prom_uzuz', ...
   'prom_uxuy', ...
   'prom_uxuz', ...
   'prom_uyuz'};
   
   groups = fields(itemsByGroup);
   
   items2 = {};
   for i = 1:length(items)
      fldnm = items{i};
      if ~ isempty(intersect(groups,{fldnm}))
         items2 = [items2 itemsByGroup.(fldnm)];
      else
         items2 = [items2 fldnm];
      end
   end

   items = [items2 {'time', 'ifrlist'}];
   
   fnm = fullfile(self.path,self.basenm,sprintf(self.mergeformat,self.basenm));
   
   fr = io.readh5file(fnm,items);
   fldlist = fields(fr);
   for i = 1:length(fldlist)
     fld = fldlist{i};
     if isempty(find(strcmp(properties(self),fld)))
        p = addprop(self,fld);
     end
     self.(fld) = fr.(fld);
   end
   
   return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function append(self, other)

   notappend = {'basenm','path'};
   flds1 = fields(self);
   flds2 = fields(self);
   flds = intersect(flds1,flds2);
   flds = setdiff(flds,notappend);
   ifr1  = self.ifrlist;
   ifr2  = other.ifrlist;
   ifr = [ ifr1 ifr2];
   self.ifrlist = ifr;
   
   for i = 1:length(flds)
      fld = flds{i};
      self.(fld) = [self.(fld) other.(fld)];
   end
   
   return
end

function complementmerged(self)

   notappend = {'basenm','path','ifrlist'};
   flds = fields(self);
   flds = setdiff(flds,notappend);

   ifr1  = self.ifrlist;
   ifr2  = self.getifr();
   ifr2 = setdiff(ifr2, ifr1);;

   ifrlist = [ ifr1 ifr2];
   
   notifr = [];
   for ifr = ifr2
      fnm = sprintf(self.fileformat, self.basenm, ifr);
      fnm = fullfile(self.path, self.basenm, fnm);
      try
         [h5f,failitems] = io.readh5file(fnm,flds);
         for i = 1:length(flds)
            fld = flds{i};
            self.(fld) = [self.(fld) h5f.(fld)'];
         end
      catch
         notifr = [notifr ifr];
      end
   end

   self.ifrlist = setdiff(ifrlist, notifr);
   
   return
end


end
end
