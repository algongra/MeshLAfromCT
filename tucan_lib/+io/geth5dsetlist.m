function geth5dsetlist(fnm,varargin)

%sizestruct = struct;
%misc.assigndefaults(varargin{:});

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

ii = find(strcmp(varargin, '-sizestruct'));
if ~isempty(ii)
   sizestruct = varargin{ii+1};
   iilist = setdiff(iilist, [ii ii+1]);
else
   sizestruct = struct;
end

for i0 = iilist
   varnm = varargin{i0};
   if ~isempty(find(strcmp(fields(sizestruct), varnm)))
      assignin('caller',[prefix varnm suffix],io.geth5dset(fnm,varnm, 'sizestruct', sizestruct.(varnm))); 
   else
      assignin('caller',[prefix varnm suffix],io.geth5dset(fnm,varnm)); 
   end
end

return
end
