function itemsByGroup = getgroups(self, ndim)

coords = {'x','y','z'};

% buffers
buffers = {'p'};
for i = 1:ndim
   myvar = ['u' coords{i}];
   buffers = [buffers myvar]; 
end
if self.nps > 0
   buffers = [buffers 'ps'];
end

% mesh
mesh = {}; 
for i = 1:ndim
   myvar1 = ['n'  coords{i}];
   myvar3 = [coords{i} '0'];
   myvar4 = [coords{i} 'f'];
   myvar = [coords{i} 'p'];
   mesh = [mesh myvar]; 
   for j = 1:ndim
      myvar = [coords{i} 'u' coords{j}];
      mesh = [mesh myvar]; 
   end
   myvar = ['d' coords{i} 'r'];
   mesh = [mesh myvar myvar3 myvar4 myvar1]; 
end

% bc
bc = {};
myvar = 'p';
for idir = 1:ndim
   for ipos = 1:2
   myvar2 = ['bc' myvar num2str(idir) num2str(ipos)];
   bc = [bc myvar2];
   end
end
for ivar = 1:ndim
   myvar = ['u' coords{ivar}];
   for idir = 1:ndim
      for ipos = 1:2
      myvar2 = ['bc' myvar num2str(idir) num2str(ipos)];
      bc = [bc myvar2];
      end
   end
end
if self.nps > 0
   myvar = 'ps';
   for idir = 1:ndim
      for ipos = 1:2
      myvar2 = ['bc' myvar num2str(idir) num2str(ipos)];
      bc = [bc myvar2];
      end
   end
end

% params
params = {};
for i = 1:ndim
    params = [params myvar1];
end

itemsByGroup.mesh    = mesh;
itemsByGroup.buffers = buffers;
itemsByGroup.bc      = bc ;


end
