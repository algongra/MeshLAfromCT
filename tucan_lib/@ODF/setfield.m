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
   % fr = ODF(basename);
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
