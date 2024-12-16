function [ifrlist stat] = getifrlist(self)
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
   stat = 0;
   dirnm = fullfile(self.datapath,self.basenm);
   if ~isdir(dirnm)
      stat = 1;
      ifrlist = [];
   else
      [ifrlist stat] = misc.getifrlist(dirnm, self.nf.regex_ifr, self.nf.ifrpos);
   end
   self.ifrlist = ifrlist;
end
