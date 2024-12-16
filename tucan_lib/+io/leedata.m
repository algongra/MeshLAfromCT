function mystruct = leedata(fname,varargin); 
%function mystruct = leefield(fname); 
%function mystruct = leefield(fname,varlooked); 
%
%  example:
%     fr = leedata('../data/re100_dt0.010_nx512.001',{'ux','uy','p'})
%
%     routine to save simulation data
%     Script to read simulation data
%
%     Allows the user to read data saved with te routine savefile.f
%      without modify this script 
%
%     an example of the file savefile is shown
%
%     EXAMPLE
%     -------
%      character(500) varlist 
%
%      write(varlist,*) 'time,d*1 ; Re,d*1; nx,i*1; ny,i*1;'//
%     .'ux,d*nx*ny'
%
%      open(110,file=filename,form='unformatted')
%      !
%      write(110) 500,';',',','*'
%      !
%      write(110) varlist
%      !
%      write(110) time, Re, nx, ny, ux
%      close(110)
%
%
varlooked='*';
found=0;
defaultvarargin={varlooked};
for i=1:length(varargin)
  if ~isempty(varargin{i})
    defaultvarargin{i} = varargin{i};
  end
end
varlooked=defaultvarargin{1};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
check=true;
fid = fopen(fname,'r','l'); 
%
% LECTURA DE PARAMETROS PARA SEPARAR EL NOMBRE DE VARIABLES Y DEFINIR
%  SUS TAMAÑOS
%
kk=fread(fid,1,'int');
listlength=fread(fid,1,'int');
itemseparator=char(fread(fid,1,'char*1'));
nameseparator=char(fread(fid,1,'char*1'));
sizeseparator=char(fread(fid,1,'char*1'));
kk=fread(fid,1,'int');
%
%  LECTURA DE VARLIST Y VARTYPE
%
kk=fread(fid,1,'int');
varlist=char(fread(fid,[1 listlength],'char*1'));
kk=fread(fid,1,'int');
%
%  TRANSFORMACION DE VARLIST Y VARTYPE
%
varlist=strrep(varlist,' ','');
varlist=regexp(varlist,itemseparator,'split');
%
kk=fread(fid,1,'int');
%
count=0;
mystruct=struct;
for ii=1:length(varlist)
   vardesc=varlist{ii};
   vardesclist=regexp(vardesc,nameseparator,'split');
   varname=vardesclist{1};
   vartypelist=regexp(vardesclist{2},sizeseparator,'split');
   prec=char(vartypelist(1));
   sizevector=[];
   for jj=2:length(vartypelist)
      mysize=int32(str2num(char(vartypelist(jj))));
      if isempty(mysize), mysize=mystruct.(char(vartypelist(jj)));,end
      sizevector=[sizevector mysize];
   end  
   if max(strcmp(varname,varlooked))==1 | strcmp('*',varlooked) | ...
      ( prod(double(sizevector))==1 && prec=='i')
      if prec == 'd' | prec == 'REAL'
          matlabformat='double';
          %value=fread(fid,sizevector,matlabformat);
          value=fread(fid,prod(double(sizevector)),matlabformat);
   if length(sizevector) == 1
         sizevector = [1 sizevector];
   end
          value = reshape(value,sizevector);
      elseif prec == 'i' | prec == 'INTEGER'
          matlabformat='int';
          %value=fread(fid,sizevector,matlabformat);
          value=fread(fid,prod(double(sizevector)),matlabformat);
   if length(sizevector) == 1
         sizevector = [1 sizevector];
   end
          value = reshape(value,sizevector);
      elseif prec == 'l' | prec == 'INTEGER(8)'
          matlabformat='int64';
          %value=fread(fid,sizevector,matlabformat);
          value=fread(fid,prod(double(sizevector)),matlabformat);
   if length(sizevector) == 1
         sizevector = [1 sizevector];
   end
          value = reshape(value,sizevector);
      elseif prec == 'c' | prec == 'CHARACTER'
          matlabformat='char*1';
          %value=strrep(char(fread(fid,[1 sizevector],matlabformat)),' ','');
          value=char(fread(fid,[1 sizevector],matlabformat));
      end
      %
      mystruct.(varname) = value;
      if max(strcmp(varlooked,varname))==1
         found=found+1;
         if found == length(varlooked);
           check=false;
           break 
         end
      end
   else
      count=prod(double(sizevector));
      if prec == 'd' | prec == 'REAL'
          fseek(fid,count*8,0);
      elseif prec == 'i' | prec == 'INTEGER'
          fseek(fid,count*4,0);
      elseif prec == 'l' | prec == 'INTEGER(8)'
          fseek(fid,count*8,0);
      elseif prec == 'c' | prec == 'CHARACTER'
          fseek(fid,count*1,0);
      end
   end
end
%
% CHECK IF THE FILE HAS BEEN CORRECTLY READ
%
if check
    kk1=fread(fid,1,'int');
    kk2=fread(fid,1,'int');
    if length(kk1) ==0
       display(sprintf('FILE %s ENDED UNEXPECTEDLY',fname))
       pause
    end
    if length(kk2) ~=0
       display(sprintf('FILE %s NOT ENDED',fname))
       pause
    end
end
fclose(fid); 
return
end
