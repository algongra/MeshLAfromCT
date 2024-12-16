function mystr = tab2ascii2(table,wlist,varargin)

cap = [];
misc.assigndefaults(varargin{:});

% get max width per column
[nrows,ncols]=size(table);

tablecell = cell(nrows, ncols);

for i = 1:nrows
   for j = 1:ncols
 
      tablecell{i,j} = misc.text2lines(table{i,j},wlist(j));
      nr(i,j) = length(tablecell{i,j});

   end
end

i2 = 1;
i3 = 1;
for i = 1:nrows
   nn = max(nr(i,:));
   i2 = i3;
   for j = 1:ncols
      for i0 = 1:nn 

         if length(tablecell{i,j}) >= i0
            tab{i2-(j-1)*nn,j} = tablecell{i,j}{i0};
         else
            tab{i2-(j-1)*nn,j} = ' ';
         end
         i2 = i2 + 1;
      end
   end
   i3 = i3 + nn;
end

mystr = mytab.tab2ascii(tab, varargin{:});



return 
end

%_______________________________________________________________________________

function mystr = hcap(mystr,cap,ncols,nmax)


nn = sum(nmax(:));
myfrmt = ['%' num2str(nn) 's'];
mystr = [mystr sprintf(myfrmt, cap) sprintf('\n')];

for j = 1:ncols
   mystr = [mystr '~'];
   for j2 = 1:nmax(j)
      mystr = [mystr '~'];
   end
end
mystr = [mystr '~' sprintf('\n')];
return 
end
%________________________________________________________________________________

function mystr = hline(mystr,ncols,nmax)

for j = 1:ncols
   mystr = [mystr '+'];
   for j2 = 1:nmax(j)
      mystr = [mystr '-'];
   end
end
mystr = [mystr '+' sprintf('\n')];
return 
end

%________________________________________________________________________________
function mystr = row(mystr,i,table,ncols,nmax)
for j = 1:ncols
   mystr=[mystr '|'];
   myfrmt = ['%' num2str(nmax(j)-1) 's'];
   uu = table{i,j};
   mystr  = [mystr sprintf(myfrmt , uu) ' '];
end
mystr = [mystr '|' sprintf('\n')];
return
end
