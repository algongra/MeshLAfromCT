function mystruct = leefulldata(fname,varargin); 

varlooked='*';
found=0;
defaultvarargin={};
defaultvarargin{1}=varlooked;
for i=1:length(varargin)
   if ~isempty(varargin{i})
    defaultvarargin{i} = varargin{i};
   end
end
varlooked=defaultvarargin{1};

iforbeg=0;
jforbeg=0;
imatbeg=1;
jmatbeg=1;
%
samevars = {'time','timearray','dt','dx','dy','Re',...
    'nproc','ndx','ndy','nx','ny',...
    'iuxbeg_glob','iuxend_glob',...
    'juxbeg_glob','juxend_glob',...
    'iuybeg_glob','iuyend_glob',...
    'juybeg_glob','juyend_glob',...
    'ipbeg_glob','ipend_glob',...
    'jpbeg_glob','jpend_glob',...
    'iuxsize_glob',...
    'juxsize_glob',...
    'iuysize_glob',...
    'juysize_glob',...
    'ipsize_glob',...
    'jpsize_glob',...
    'basename','tol','timeframemax','iframe',...
    'CFLarray','svnrevision','source','shellscript',...
    'solverflagux','solverflaguy','solverflagp'...
    'deltapoint','lsize','testcase','x0','xf','y0','yf',...
    'norminf',...
    'xib','yib','xib0','yib0','fx','fy','lux','luy','ncount'};
%
fullvars = {'ux','uy','p'};
%
fullvars2 = {'phi','rhsp','rhsux','rhsuy','uxwk','uywk'};
ifullvars2 = {'p','p','ux','uy','ux','uy'};
%
horvars = {'xux','xuy','xp'};
ihorvars = {'ux','uy','p'};
%
vervars = {'yux','yuy','yp'};
ivervars = {'ux','uy','p'};

vervars2 = {'uxeast','uyeast'};
ivervars2 = {'ux','uy'};
%
seqvars = {'timingvalues'};
% fname may be whichever of the saved files (by any processor)
auxlist=regexp(fname,'\.','split');
fnmnocpu='';
for i = 1:length(auxlist)-1
   fnmnocpu = [fnmnocpu auxlist{i} '.'];
end

fr = io.leedata(fname,{'nproc','ndx','ndy'});
nproc = fr.nproc;ndx=fr.ndx;ndy=fr.ndy;
fullbuffcell = cell(1,length(fullvars));
fullbuffcell2= cell(1,length(fullvars2));
horbuffcell  = cell(1,length(horvars));
verbuffcell  = cell(1,length(vervars));
verbuffcell2 = cell(1,length(vervars2));
seqbuffcell  = cell(1,length(seqvars));
for i = 1:ndx
  for j = 1:ndy
    ii = j-1+(i-1)*ndy;
    fr = io.leedata([fnmnocpu sprintf('%03d',ii)],varlooked);
    % full
    for ivar = 1:length(fullvars)
        var = fullvars{ivar};
        if isfield(fr,var)
           ibeg = getfield(fr,['i' var 'beg'])+imatbeg-iforbeg;
           iend = getfield(fr,['i' var 'end'])+imatbeg-iforbeg;
           jbeg = getfield(fr,['j' var 'beg'])+jmatbeg-jforbeg;
           jend = getfield(fr,['j' var 'end'])+jmatbeg-jforbeg;
           %
           ibeg0 = getfield(fr,['i' var 'beg0'])+imatbeg-iforbeg;
           jbeg0 = getfield(fr,['j' var 'beg0'])+jmatbeg-jforbeg;
           ioff=ibeg-ibeg0;joff=jbeg-jbeg0;
           %
           buff = fullbuffcell{ivar};
           buff (ibeg:iend,jbeg:jend) = ...
            fr.(var)(1+ioff:iend-ibeg+1+ioff,1+joff:jend-jbeg+1+joff);
           fullbuffcell{ivar} = buff;
        end
    end   
    for ivar = 1:length(fullvars2)
        var = fullvars2{ivar};
        if isfield(fr,var)
           ibeg = getfield(fr,['i' ifullvars2{ivar} 'beg'])+imatbeg-iforbeg;
           iend = getfield(fr,['i' ifullvars2{ivar} 'end'])+imatbeg-iforbeg;
           jbeg = getfield(fr,['j' ifullvars2{ivar} 'beg'])+jmatbeg-jforbeg;
           jend = getfield(fr,['j' ifullvars2{ivar} 'end'])+jmatbeg-jforbeg;
           %
           ibeg0=getfield(fr,['i' ifullvars2{ivar} 'beg0'])+imatbeg-iforbeg;
           jbeg0=getfield(fr,['j' ifullvars2{ivar} 'beg0'])+jmatbeg-jforbeg;
           ioff=ibeg-ibeg0;joff=jbeg-jbeg0;
           %
           buff = fullbuffcell2{ivar};
           buff (ibeg:iend,jbeg:jend) = ...
            fr.(var)(1+ioff:iend-ibeg+1+ioff,1+joff:jend-jbeg+1+joff);
           fullbuffcell2{ivar} = buff;
        end
    end   
    % horizontal
    for ivar = 1:length(horvars)
        var = horvars{ivar};
        if isfield(fr,var)
           ibeg = getfield(fr,['i' ihorvars{ivar} 'beg'])+imatbeg-iforbeg;
           iend = getfield(fr,['i' ihorvars{ivar} 'end'])+imatbeg-iforbeg;
           ibeg0 = getfield(fr,['i' ihorvars{ivar} 'beg0'])+imatbeg-iforbeg;
           ioff=ibeg-ibeg0;
           buff = horbuffcell{ivar};
           buff (ibeg:iend) = fr.(var)(1+ioff:iend-ibeg+1+ioff);
           horbuffcell{ivar} = buff;
        end
    end   
    % vertical
    for ivar = 1:length(vervars)
        var = vervars{ivar};
        if isfield(fr,var)
           jbeg = getfield(fr,['j' ivervars{ivar} 'beg'])+imatbeg-iforbeg;
           jend = getfield(fr,['j' ivervars{ivar} 'end'])+imatbeg-iforbeg;
           jbeg0 = getfield(fr,['j' ivervars{ivar} 'beg0'])+jmatbeg-jforbeg;
           joff=jbeg-jbeg0;
           buff = verbuffcell{ivar};
           buff (jbeg:jend) = fr.(var)(1+joff:jend-jbeg+1+joff);
           verbuffcell{ivar} = buff;
        end
    end   
    % vertical
    for ivar = 1:length(vervars2)
        var = vervars2{ivar};
        if isfield(fr,var)
           jbeg = getfield(fr,['j' ivervars2{ivar} 'beg'])+imatbeg-iforbeg;
           jend = getfield(fr,['j' ivervars2{ivar} 'end'])+imatbeg-iforbeg;
           buff = verbuffcell2{ivar};
           buff (jbeg:jend) = fr.(var)(1:jend-jbeg+1);
           verbuffcell2{ivar} = buff;
        end
    end   
    % seq
    for ivar = 1:length(seqvars)
        var = seqvars{ivar};
        if isfield(fr,var)
           buff = seqbuffcell{ivar};
           buff{ii+1} = fr.(var);
           seqbuffcell{ivar} = buff;
        end
    end   
  end
end
mystruct = struct;

% full
for ivar = 1:length(fullvars);
   var = fullvars{ivar};
   buff = fullbuffcell{ivar};
   mystruct.(var)=buff;
end
% full
for ivar = 1:length(fullvars2);
   var = fullvars2{ivar};
   buff = fullbuffcell2{ivar};
   mystruct.(var)=buff;
end
% hor
for ivar = 1:length(horvars);
   var = horvars{ivar};
   buff = horbuffcell{ivar};
   mystruct.(var)=buff;
end
% ver
for ivar = 1:length(vervars);
   var = vervars{ivar};
   buff = verbuffcell{ivar};
   mystruct.(var)=buff;
end
% ver
for ivar = 1:length(vervars2);
   var = vervars2{ivar};
   buff = verbuffcell2{ivar};
   mystruct.(var)=buff;
end
% full
for ivar = 1:length(samevars);
   var = samevars{ivar};
   if isfield(fr,var)
      mystruct.(var)=fr.(var);
   end
end
% seq
for ivar = 1:length(seqvars);
   var = seqvars{ivar};
   buff = seqbuffcell{ivar};
   mystruct.(var)=buff;
end
return
end
