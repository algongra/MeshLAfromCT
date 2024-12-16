function bytes = getdiskspace(self,varargin)
   % @author M.Moriche
   % @date 21-10-2014
   % @brief Method to get the disk space in bytes used by the simulation
   %
   % @details
   %
   % Makes use of UNIX built-in application du
   %
   % @code
   % fr = Frame2D(basenm, 'path', datapath);
   % bb = fr.getdiskspace();   % bytes
   % Gb = bb/(2^30);           % Gyga bytes
   % @endcode
  
   mypatt = '(.*\n)*(\d*)(.*)total';

   bytes = 0;
   [st,rs] = system(sprintf('du -sb --total %s/%s/*%s', self.datapath, self.basenm,self.nf.ext));

   if st == 0
      a = regexp(rs, mypatt, 'tokens');
      bb = str2num(a{1}{2});
      bytes = bytes + bb;
   else
      display('error... ')
      display(rs)
   end
end
