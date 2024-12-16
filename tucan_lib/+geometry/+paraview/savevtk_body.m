function savevtk(ver,fac,filename); 

fid=fopen(filename,'wt'); 

fprintf(fid,'%s\n','# vtk DataFile Version 2.0'); 
fprintf(fid,'%s\n','Send complaints about this file to guillem'); 
%---
fprintf(fid,'%s\n','ASCII'); 
%fprintf(fid,'%s\n','BINARY'); 
%---
fprintf(fid,'%s\n','DATASET POLYDATA'); 
fprintf(fid,'%s %i %s\n','POINTS',size(ver,1), 'float'); 

for i=1:size(ver,1); 
    fprintf(fid,'%f %f %f\n',ver(i,1),ver(i,2),ver(i,3)); 
end
fprintf(fid,'%s %i %i\n','POLYGONS',size(fac,1), size(fac,1)*(size(fac,2)+1)); 
for i=1:size(fac,1); 
    %fprintf(fid,'%i %i %i %i %i\n',4,fac(i,1)-1,fac(i,2)-1,fac(i,3)-1,fac(i,4)-1); 
    fprintf(fid,'%i ',size(fac,2)); 
    for j = 1:size(fac,2)
       fprintf(fid,'%i ',fac(i,j)-1);
    end
    fprintf(fid,'\n');
end
fclose(fid); 
return
