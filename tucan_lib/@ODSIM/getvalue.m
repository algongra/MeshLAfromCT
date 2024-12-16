function value = getvalue(self, varnm)


%misc.assigndefaults(varargin{:});
ifrlist = self.ifrlist;

% get dimensions
nnlist = zeros(length(ifrlist), 1);
n0list = zeros(length(ifrlist), 1);
for i0 = 1:length(ifrlist)
   ifr = ifrlist(i0);
   self.nf.ifr = ifr;
   fnm = self.nf.getfullfilename();
   dims = io.getdimensions(fnm, varnm);
   nnlist(i0) = prod(dims);
   n0list(i0) = dims(1);
   n2list(i0) = prod(dims(2:end));
end

n0 = sum(n0list);
nn = sum(nnlist);
value = zeros(nn,1);

ioff = length(ifrlist);

ia0 = 1;
for i0 = 1:length(ifrlist)

   ifr = ifrlist(i0);
   self.nf.ifr = ifr;
   fnm = self.nf.getfullfilename();

   vv = io.geth5dset(fnm, varnm);
   vv2 = reshape(vv, n0list(i0), n2list(i0));

   ia = ia0;
   for i2 = 1:n2list(i0)
      ib = ia + n0list(i0)-1;
      value(ia:ib) = vv2(:,i2);
      ia = ia + n0;
   end
   ia0 = ia0 + n0list(i0);

end


if length(dims) == 1
   newdims = [n0 1];
else
   newdims = [n0 dims(2:end)]
end

value = reshape(value, newdims);

return
end




