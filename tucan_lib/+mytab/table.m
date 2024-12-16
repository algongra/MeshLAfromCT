function tab = table(basenmlist, paramlist, varargin)

asciiflag = false;
ifr = 'last';
path = '.';
misc.assigndefaults(varargin{:});

fr = Frame2D(basenmlist{1}, varargin{:});
fnm = fr.getfullfilename();
io = io.readh5file(fnm, {'ndim'});

try
   if io.ndim == 2
      mytab.data;
   elseif io.ndim == 3
      mytab.data3D;
   end
catch
   mytab.data;
end
if asciiflag
   head = 'headascii';
   patt = 'patternascii';
   editbasenm = @(basenm) upper(basenm);
else
   head = 'head';
   patt = 'pattern';
   editbasenm = @(basenm) ['{\verb ' upper(basenm) ' }'];
end


% GET FIELD LIST ----------------------------------------------
fldlist = {};
for i = 1:length(paramlist)
   param = paramlist{i};
   if isfield(auxByHeader, param)
      for j = 1:length(auxByHeader.(param).fields)
         item = auxByHeader.(param).fields{j};
         a = regexp(item, '(\w*)\((.*)\)','tokens');
         if length(a) == 0
            fldlist = [fldlist item];
         else
            fldlist = [fldlist a{1}(1:2:end)'];
         end
      end
   else
      fldlist = [fldlist param];
   end
end

% HEADER
tab = cell(1,length(paramlist) + 1);
row = {'case'};
for i = 1:length(paramlist)
   param = paramlist{i};
   if isfield(auxByHeader, param)
      row = [row auxByHeader.(param).(head)];
   else
      row = [row param];
   end
end

tab(1,:) = row(:);

% READ BASENM'S and fill rows

for k = 1:length(basenmlist)
   basenm = basenmlist{k};
   fr = Frame2D(basenm, 'path', path, 'ifr', ifr);
   fr.readfields('items',fldlist(:));
   row = {editbasenm(basenm)};
   for i = 1:length(paramlist)
      param = paramlist{i};

      mylist = {};
      if isfield(auxByHeader, param)
         for j = 1:length(auxByHeader.(param).fields)
            item = auxByHeader.(param).fields{j};
            a = regexp(item, '(\w*)\((.*)\)','tokens');
            if length(a) == 0
               mylist = [mylist getfield(fr,item)];
            else
               uu = a{1}(:);
               for ii = 1:2:length(uu)
                  uu2 = getfield(fr,uu{ii});
                  mylist = [mylist uu2(str2num(uu{ii+1}))];
               end
            end
         end
         if ismember('operation',fields(auxByHeader.(param)))
            eval(auxByHeader.(param).operation);
            row = [row sprintf(auxByHeader.(param).(patt), mylist{:})];
         else
            row = [row sprintf(auxByHeader.(param).(patt), mylist{:})];
         end
      else
         item = getfield(fr,param);
            if ischar(item)
               row = [row item];
            elseif isnumeric(item)
               row = [row num2str(item)];
            end
      end
   end

   tab(k+1,:) = row(:);
end


return
end
