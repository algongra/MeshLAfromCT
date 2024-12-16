function [ifr varargout]= getifr(self, ifr)

[ifrlist stat] = self.getifrlist();


if stat == 0

   if isstr(ifr)
      if strcmp(ifr,'last')
         ifr = ifrlist(end);
      elseif strcmp(ifr,'first')
         ifr = ifrlist(1);
      end
   else
      ii = find(ifrlist == ifr);
      if isempty(ii)
         misc.warningmessage('IFR not available')
         stat = 1;
      else
         ifr = ifrlist(ii);
      end
   end

end

varargout{1} = stat;

return
end

