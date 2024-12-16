function WriteToVTR(fname,x,y,z,varargin)
% Writes a VTR rectangular grid file to fname (data are written in binary
% form using Float32 type
%   WriteToVTR(fname,x,y,z,varargin)
%   fname [string]: file name
%   x,y,z [1D arrays]: the length of these arrays must match the fields
%       input later
%   the fields must be input with a string defining name and type of field
%       and the 3D array itself, size(f) = [nx,ny,nz]+
%
%   scalar integer field must end with a 'i' character, e.g.:
%       WriteToVTR('fname.vtr',x,y,z,'field1_i',field1)
%
%   vector real field must end with a 'v' character, and the three
%   directions must be input as separate arrays one after the other one.
%       WriteToVTR('fname.vtr',x,y,z,'field2_v',field2x,field2y,field2z)
%
%   scalar real field must NOT end with a 'i' or 'v' character, e.g.:
%       WriteToVTR('fname.vtr',x,y,z,'field3',field3)
%   
nx = length(x);
ny = length(y);
nz = length(z);
fid = fopen(fname,'w');%ASCII file header
offset = 0;
fprintf(fid, '<?xml version="1.0"?>\n<VTKFile type="RectilinearGrid" version="0.1" byte_order="LittleEndian">\n');
fprintf(fid, ' <RectilinearGrid WholeExtent=" %d %d %d %d %d %d">\n',0,nx-1,0,ny-1,0,nz-1);
fprintf(fid, '  <Piece Extent=" %d %d %d %d %d %d">\n',0,nx-1,0,ny-1,0,nz-1);
% Coordinates
fprintf(fid, '   <Coordinates>\n');
fprintf(fid, '    <DataArray type="Float32" Name="X" format="appended" offset="%i">',offset);
offset = offset+(nx+1)*4;
fprintf(fid, '\n    </DataArray>\n');
fprintf(fid, '    <DataArray type="Float32" Name="Y" format="appended" offset="%i">',offset);
offset = offset+(ny+1)*4;
fprintf(fid, '\n    </DataArray>\n');
fprintf(fid, '    <DataArray type="Float32" Name="Z" format="appended" offset="%i">',offset);
offset = offset+(nz+1)*4;
fprintf(fid, '\n    </DataArray>\n');
fprintf(fid, '   </Coordinates>\n');

% Data values
fprintf(fid, '   <PointData>\n');
ind = 1;
while ind < length(varargin)
    select_case = lower(varargin{ind});
    % Use last two character if possible, otherwise treat as real scalar
    if length(select_case)>2
        select_case = select_case(end-1:end);
    end
    switch select_case
        case '_v'
            fprintf(fid,'    <DataArray type="Float32" Name="%s" NumberOfComponents="3" format="appended" offset="%i">',varargin{ind},offset);
            fprintf(fid, '\n    </DataArray>\n');
            offset = offset+(3*nx*ny*nz)*4+4;
            ind = ind + 4;
        case '_i'
            fprintf(fid,'    <DataArray type="Int16" Name="%s" NumberOfComponents="1" format="appended" offset="%i">',varargin{ind},offset);
            offset = offset+(nx*ny*nz)*2+4;
            fprintf(fid, '\n    </DataArray>\n');
            ind = ind + 2;            
        otherwise
            fprintf(fid,'    <DataArray type="Float32" Name="%s" NumberOfComponents="1" format="appended" offset="%i">',varargin{ind},offset);
            offset = offset+(nx*ny*nz)*4+4;
            fprintf(fid, '\n    </DataArray>\n');
            ind = ind + 2;            
    end
end
fprintf(fid, '   </PointData>\n');
fprintf(fid, '  </Piece>\n');
fprintf(fid, ' </RectilinearGrid>\n');

% All raw data
fprintf(fid,'  <AppendedData encoding="raw">\n_');
% Coordinate Blocks
fwrite(fid,nx*4,'uint32');
fwrite(fid,x,'float32');
fwrite(fid,ny*4,'uint32');
fwrite(fid,y,'float32');
fwrite(fid,nz*4,'uint32');
fwrite(fid,z,'float32');

ind = 1;
while ind < length(varargin)
    select_case = lower(varargin{ind});
    select_case = select_case(end);
    switch select_case
        case 'v'
            VelWrite = zeros(3,nx,ny,nz);
            VelWrite(1,:,:,:) = reshape(varargin{ind+1},[1 nx ny nz]);
            VelWrite(2,:,:,:) = reshape(varargin{ind+2},[1 nx ny nz]);
            VelWrite(3,:,:,:) = reshape(varargin{ind+3},[1 nx ny nz]);
            fwrite(fid,3*nx*ny*nz*4,'uint32');
            fwrite(fid,VelWrite,'float32');
            ind = ind + 4;
        case 'i'
            fwrite(fid,nx*ny*nz*2,'uint32');
            fwrite(fid,varargin{ind+1},'int16');
            ind = ind + 2;
        otherwise
            fwrite(fid,nx*ny*nz*4,'uint32');
            fwrite(fid,varargin{ind+1},'float32');
            ind = ind + 2;
    end
end
fprintf(fid,'\n  </AppendedData>\n');
fprintf(fid, '</VTKFile>');
fclose(fid);
