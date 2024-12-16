function readfields(self,items, overrides)
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
   % fr.readfields({'xib','yib'});
   % fr.readfields({'ux','uy'},false);
   % fr.readfields({'params','mesheu'},true);
   % @endcode
   %

   % flag to override existing items
   fnm = sprintf(self.frmt, self.basenm,self.ifr,self.ext);
   fnm = fullfile(self.path, self.basenm, fnm);
   [fr,failitems] = io.readh5file(fnm,items,'sizebyvar',self.sizebyvar);
   if length(failitems) > 0
      %display('    items not found:');
      failitems;
   end

   fldlist = fields(fr);
   for i =1:length(fldlist)
      fld = fldlist{i};
      if ~isempty(getfield(fr,fld))
         self.setfield(fld,getfield(fr,fld),'overrides',overrides);
      end
   end

end 
