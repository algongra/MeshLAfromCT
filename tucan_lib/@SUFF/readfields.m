function readfields(self, fnm, items)

groups = fields(self.itemsByGroup);
items2 = {};
for i = 1:length(items)
   fldnm = items{i};
   if ~ isempty(intersect(groups,{fldnm}))
      items2 = [items2 self.itemsByGroup.(fldnm)];
   else
      items2 = [items2 fldnm];
   end
end
items = items2;

slabvars = intersect(items, fields(self.sizebyvar));

restvars = setdiff(items, slabvars);

% READ SLAB VARS

for i0 = 1:length(slabvars)
   varnm = slabvars{i0};
   sizestruct = self.sizebyvar.(varnm);
   self.(varnm) = io.geth5dset(fnm, varnm, 'sizestruct', self.sizebyvar.(varnm));
end

% READ NON-SLAB VARS

for i0 = 1:length(restvars)
   varnm = restvars{i0};
   self.(varnm) = io.geth5dset(fnm, varnm);
end

return
end
