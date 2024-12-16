function getlist(self, varargin)

iilist = 1:length(varargin);
ii = find(strcmp(varargin, '-prefix'));
if ~isempty(ii)
   prefix = varargin{ii+1};
   iilist = setdiff(iilist, [ii ii+1]);
else
   prefix = '';
end

ii = find(strcmp(varargin, '-suffix'));
if ~isempty(ii)
   suffix = varargin{ii+1};
   iilist = setdiff(iilist, [ii ii+1]);
else
   suffix = '';
end

for i0 = iilist
   varnm = varargin{i0};
   assignin('caller',[prefix varnm suffix],self.getvalue(varnm)); 
end

return
end

