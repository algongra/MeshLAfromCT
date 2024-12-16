function bytes = getdiskspace(self)
   % @author M.Moriche
   % @date 21-10-2014
   % @brief Method to get the disk space in bytes used by the simulation
   %
   % @details
   %
   % @code
   % fr = ODF2D(basenm, 'path', path);
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

