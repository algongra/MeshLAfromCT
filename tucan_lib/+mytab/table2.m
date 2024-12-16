function tab = table2(basenm, paramlist, varargin)

asciiflag = false;
ifr = 'last';
path = '.';
misc.assigndefaults(varargin{:});

fr = Frame2D(basenm, varargin{:});
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
   editbasenm = @(basenm) upper(strrep(basenm,'_','\_'));
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

fr = Frame2D(basenm, 'path', path, 'ifr', ifr);
fr.readfields('items',fldlist);

% HEADER
tab = cell(1,2);
tab(1,:) = {'Parameter','Value'};
for i = 1:length(paramlist)
   param = paramlist{i};
   if isfield(auxByHeader, param)
      mylist = {};
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
         row = sprintf(auxByHeader.(param).(patt), mylist{:});
      else
         row = sprintf(auxByHeader.(param).(patt), mylist{:});
      end

      %tab(i+1,:) = {strrep(param,'_',' '), sprintf(auxByHeader.(param).(patt), mylist{:})};
      tab(i+1,:) = {['{' auxByHeader.(param).head '}'], sprintf(auxByHeader.(param).(patt), mylist{:})};
   else
      item = getfield(fr,param);
      if ischar(item)
         row = {strrep(param,'_',' '),  item};
      elseif isnumeric(item)
         row = {strrep(param,'_',' '), num2str(item)};
      end
      tab(i+1,:) = row;
   end
end

return
end
