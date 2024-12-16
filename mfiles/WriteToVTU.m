function WriteToVTU(fname,pts,cells,varargin)
% Write unstructured mesh with data to .vtu xml file
% Usage:
%   WriteToVTU(fname,pts,cells,varargin)
%       fname [string]  : file name
%       pts   [3 x np]  : location of vertex points in 3D
%       cells [nv x nc] : connectivity matrix for each face (for a triangle
%                       nv would be 3)
%       varargin        : cell values, first a string with the name to use
%                       to save the variable, followed by the array itself,
%                       e.g: 'data',[1 0.2 3.5 4 5.2] by default float32
%                       the last two character of the string can be used to
%                       specify a vector ('_v') or an integer quantify
%                       ('_i'). e.g.: 'vel_v',[1 2 3],[0 2 3],[5 6 0]
%
% Created by Lorenzo Rossini <lorenz.ross@gmail.com>
% May 2016
pts_size = size(pts);
cells_size = size(cells);

% Make sure the triangle and the points are not transposed
if cells_size(1)>cells_size(2)
    cells = cells';
    nc = cells_size(1);
else
    nc = cells_size(2);
end
if pts_size(1)~=3 && pts_size(2)~=3
    error('The points array needs the 3 [x y z] components in space');
elseif pts_size(1)>pts_size(2)
    pts = pts';
    np = pts_size(1);
else
    np = pts_size(2);
end

% Check if a cell type is provided, otherwise use a triangle (5) by default
cell_type_i_arg = 0;
write_cell_data = true;
for i_arg=1:length(varargin)
    if strcmpi(varargin{i_arg},'cell_type') ||...
            strcmpi(varargin{i_arg},'cell') ||...
            strcmpi(varargin{i_arg},'celltype')
        cell_type_i_arg = i_arg;
    end
end
if cell_type_i_arg
    cell_types = varargin{cell_type_i_arg+1};
    cell_offsets = varargin{cell_type_i_arg+2};
    varargin = {varargin{1:cell_type_i_arg-1} , varargin{cell_type_i_arg+3:end}};
    if all(cell_types==1)
        write_cell_data = false;
    end
else
    cell_types = 5*ones([1 nc]);
    cell_offsets = (3:3:3*nc);
end

% Find out if each field is a cell or point data field
field_type = zeros(length(varargin),1);
for ind = 1:length(varargin)
    if ischar(varargin{ind}(:))
        field_type(ind) = 0;
    elseif length(varargin{ind}(:)) == nc
        field_type(ind) = 2;
        field_type(ind-1) = 2;
    elseif length(varargin{ind}(:)) == np
        field_type(ind) = 1;
        field_type(ind-1) = 1;
    else
        if ischar(varargin{ind-1})
            error([' Field "' varargin{ind-1} '" (n.' num2str(ind) ') has inconsistent dimensions!']);
        else
            error([' Field n.' num2str(ind) ' has inconsistent dimensions!']);
        end
    end
end

%ASCII file header
fid = fopen(fname,'w');
offset = 0;
fprintf(fid, '<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n');
fprintf(fid, ' <UnstructuredGrid>\n');
fprintf(fid, '  <Piece NumberOfPoints="%i" NumberOfCells="%i">\n',np,nc);

%% Point data header
fprintf(fid, '   <PointData>\n');
ind = 1;
while ind < length(varargin)
    select_case = lower(varargin{ind});
    % Use last two character if possible, otherwise treat as real scalar
    if length(select_case)>2
        select_case = select_case(end-1:end);
    end
    % Check the data is a point data field
    if length(varargin{ind+1})==np
        switch select_case
            case '_v'
                fprintf(fid,'    <DataArray type="Float32" Name="%s" NumberOfComponents="3" format="appended" offset="%i">',varargin{ind},offset);
                fprintf(fid, '\n    </DataArray>\n');
                offset = offset+(3*np)*4+4;
                ind = ind + 4;
            case '_i'
                fprintf(fid,'    <DataArray type="Int32" Name="%s" NumberOfComponents="1" format="appended" offset="%i">',varargin{ind},offset);
                offset = offset+(np)*4+4;
                fprintf(fid, '\n    </DataArray>\n');
                ind = ind + 2;
            otherwise
                fprintf(fid,'    <DataArray type="Float32" Name="%s" NumberOfComponents="1" format="appended" offset="%i">',varargin{ind},offset);
                offset = offset+(np)*4+4;
                fprintf(fid, '\n    </DataArray>\n');
                ind = ind + 2;
        end
    else
        ind = ind + 1;
    end
end
fprintf(fid, '   </PointData>\n');

%% Cell data header
fprintf(fid, '   <CellData>\n');
ind = 1;
while ind < length(varargin)
    select_case = lower(varargin{ind});
    % Use last two character if possible, otherwise treat as real scalar
    if length(select_case)>2
        select_case = select_case(end-1:end);
    end
    if length(varargin{ind+1})==nc && write_cell_data
        switch select_case
            case '_v'
                fprintf(fid,'    <DataArray type="Float32" Name="%s" NumberOfComponents="3" format="appended" offset="%i">',varargin{ind},offset);
                fprintf(fid, '\n    </DataArray>\n');
                offset = offset+(3*nc)*4+4;
                ind = ind + 4;
            case '_i'
                fprintf(fid,'    <DataArray type="Int32" Name="%s" NumberOfComponents="1" format="appended" offset="%i">',varargin{ind},offset);
                offset = offset+(nc)*4+4;
                fprintf(fid, '\n    </DataArray>\n');
                ind = ind + 2;
            otherwise
                fprintf(fid,'    <DataArray type="Float32" Name="%s" NumberOfComponents="1" format="appended" offset="%i">',varargin{ind},offset);
                offset = offset+(nc)*4+4;
                fprintf(fid, '\n    </DataArray>\n');
                ind = ind + 2;
        end
    else
        ind = ind + 1;        
    end
end
fprintf(fid, '   </CellData>\n');

%% Mesh Data header
fprintf(fid, '   <Points>\n');
fprintf(fid, '    <DataArray type="Float32" Name="Coordinates" NumberOfComponents="3" format="appended" offset="%i">',offset);
offset = offset+(3*np+1)*4;
fprintf(fid, '\n    </DataArray>\n');
fprintf(fid, '   </Points>\n');

fprintf(fid, '   <Cells>\n');
if nc>0
    fprintf(fid, '    <DataArray type="Int32" Name="connectivity" format="appended" offset="%i">',offset);
    offset = offset+(cell_offsets(1)*nc)*4+4;
    fprintf(fid, '\n    </DataArray>\n');
    fprintf(fid, '    <DataArray type="Int32" Name="offsets" format="appended" offset="%i">',offset);
    offset = offset+(nc)*4+4;
    fprintf(fid, '\n    </DataArray>\n');
    fprintf(fid, '    <DataArray type="UInt8" Name="types" format="appended" offset="%i">',offset);
    offset = offset+(nc)*2+4;
    fprintf(fid, '\n    </DataArray>\n');
end
fprintf(fid, '   </Cells>\n');

fprintf(fid, '  </Piece>\n');
fprintf(fid, ' </UnstructuredGrid>\n');

%%%%%%%%%%%%%%% All raw data %%%%%%%%%%%%%%%
fprintf(fid,'  <AppendedData encoding="raw">\n_');
%% Point Data first
ind = 1;
while ind < length(varargin)
    select_case = lower(varargin{ind});
    if length(select_case)>2
        select_case = select_case(end-1:end);
    end
    if length(varargin{ind+1})==np
        switch select_case
            case '_v'
                VelWrite = zeros(3,np);
                VelWrite(1,:) = reshape(varargin{ind+1},[1 np]);
                VelWrite(2,:) = reshape(varargin{ind+2},[1 np]);
                VelWrite(3,:) = reshape(varargin{ind+3},[1 np]);
                fwrite(fid,3*np*4,'uint32');
                fwrite(fid,VelWrite,'float32');
                ind = ind + 4;
            case '_i'
                fwrite(fid,np*4,'uint32');
                fwrite(fid,varargin{ind+1},'int32');
                ind = ind + 2;
            otherwise
                fwrite(fid,np*4,'uint32');
                fwrite(fid,varargin{ind+1},'float32');
                ind = ind + 2;
        end
    else
        ind = ind + 1;        
    end
end

%% Cell Data second
ind = 1;
while ind < length(varargin)
    select_case = lower(varargin{ind});
    if length(select_case)>2
        select_case = select_case(end-1:end);
    end
    if length(varargin{ind+1})==nc && write_cell_data
        switch select_case
            case '_v'
                VelWrite = zeros(3,nc);
                VelWrite(1,:) = reshape(varargin{ind+1},[1 nc]);
                VelWrite(2,:) = reshape(varargin{ind+2},[1 nc]);
                VelWrite(3,:) = reshape(varargin{ind+3},[1 nc]);
                fwrite(fid,3*nc*4,'uint32');
                fwrite(fid,VelWrite,'float32');
                ind = ind + 4;
            case '_i'
                fwrite(fid,nc*4,'uint32');
                fwrite(fid,varargin{ind+1},'int32');
                ind = ind + 2;
            otherwise
                fwrite(fid,nc*4,'uint32');
                fwrite(fid,varargin{ind+1},'float32');
                ind = ind + 2;
        end
    else
        ind = ind + 1;        
    end
end

%% Mesh Data finally
fwrite(fid,4*3*np,'uint32');
fwrite(fid,pts,'float32');
if nc>0
    fwrite(fid,4*cell_offsets(1)*nc,'uint32');
    fwrite(fid,cells-1,'int32');
    fwrite(fid,4*nc,'uint32');
    fwrite(fid,cell_offsets,'int32');
    fwrite(fid,nc,'uint32');
    fwrite(fid,cell_types,'uint8');
end
fprintf(fid,'\n  </AppendedData>\n');
fprintf(fid, '</VTKFile>');
fclose(fid);
