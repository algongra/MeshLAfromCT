function mystr = tab2ascii(table,varargin)

cap = [];
misc.assigndefaults(varargin{:});

% get max width per column
[nrows,ncols]=size(table);
nmax = zeros(ncols, 1);
for i = 1:nrows
   for j = 1:ncols
      n=length(table{i,j});
      if n > nmax(j)
         nmax(j)=n;
      end
   end
end
nmax = nmax + 2;

mystr = '';
mystr = hlinetop(mystr,ncols,nmax);
mystr = row(mystr,1,table,ncols,nmax);
mystr = hline(mystr,ncols,nmax);
for i = 2:nrows
   mystr = row(mystr,i,table,ncols,nmax);
end
mystr = hlinebot(mystr,ncols,nmax);

if ~isempty(cap)
   mystr = hcap(mystr,cap,ncols,nmax);
end

return 
end

%________________________________________________________________________________

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
   if j == 1
      %mystr = [mystr '├'];
      mystr = [mystr '|'];
   else
      %mystr = [mystr '┼'];
      mystr = [mystr '|'];
   end
   for j2 = 1:nmax(j)
      %mystr = [mystr '─'];
      mystr = [mystr '-'];
   end
end
%mystr = [mystr '|' sprintf('\n')];
mystr = [mystr '|' sprintf('\n')];
return 
end
function mystr = hlinebot(mystr,ncols,nmax)

for j = 1:ncols
   if j == 1
      %mystr = [mystr '└'];
      mystr = [mystr '+'];
   else
      %mystr = [mystr '┴'];
      mystr = [mystr '-'];
   end
   for j2 = 1:nmax(j)
      %mystr = [mystr '─'];
      mystr = [mystr '-'];
   end
end
%mystr = [mystr '┘' sprintf('\n')];
mystr = [mystr '+' sprintf('\n')];
return 
end

function mystr = hlinetop(mystr,ncols,nmax)

for j = 1:ncols
   if j == 1
      %mystr = [mystr '┌'];
      mystr = [mystr '+'];
   else
      %mystr = [mystr '┬'];
      mystr = [mystr '+'];
   end
   for j2 = 1:nmax(j)
      %mystr = [mystr '─'];
      mystr = [mystr '-'];
   end
end
%mystr = [mystr '┐' sprintf('\n')];
mystr = [mystr '+' sprintf('\n')];
return 
end

%________________________________________________________________________________
function mystr = row(mystr,i,table,ncols,nmax)
for j = 1:ncols
   %mystr=[mystr '│'];
   mystr=[mystr '|'];
   myfrmt = ['%' num2str(nmax(j)-1) 's'];
   uu = table{i,j};
   mystr  = [mystr sprintf(myfrmt , uu) ' '];
end
%mystr = [mystr '│' sprintf('\n')];
mystr = [mystr '|' sprintf('\n')];
return
end
