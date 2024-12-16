function var = geth5dset(fnm, varnm, varargin)
%%
%% @author M.Moriche                                                      %       
%% @brief Function to get dataset $varnm from a file saved with HDF5 lib  %
%% @date 18-05-2015 by M.Moriche \n                                       %
%%       Created                                                          %
%%                                                                        %
%% @details                                                               %
%%                                                                        %
%% Uses low level function in HDF5 library. High level libraries are      %
%%  not available in all Matlab releases.                                 %
%%                                                                        %
%% Mandatory arguments                                                    %
%%  - fnm: HDF5 file name to be read.                                     %
%%  - varnm: Name of the variable to read                                 %
%%                                                                        %
%% Optional arguments.                                                    %
%%  - sizestruct =[]: struct to define the slicing of the dataset         %
%%            keys = {offset, stride, block, count}                       %
%%  - blank = true: boolean to indicate if character strings should be    %
%%                  deblanked.                                            %
%% Example:                                                               %
%%                                                                        %
%% @code                                                                  %
%% kz = io.geth5dset('/data2/mmoriche/myfile.h5', 'kz');                  %
%% siz.offset = [ 0  0];                                                  %
%% siz.stride = [ 1  1];                                                  %
%% siz.block  = [ 1  1];                                                  %
%% siz.count  = [32 32];                                                  %
%% ux=io.geth5dset('/data2/mmoriche/myfile.h5','ux','sizestruct', siz);   %
%% @endcode                                                               %
%%                                                                        %
%%                                                                        %
%%                                                                        %

sizestruct.offset = [];
sizestruct.stride = [];
sizestruct.block  = [];
sizestruct.count  = [];

blank     = true;
misc.assigndefaults(varargin{:});

failitems = {};

h5fid = H5F.open(fnm,'H5F_ACC_RDONLY','H5P_DEFAULT');

did = H5D.open(h5fid,varnm);

if ~isempty(sizestruct.offset)
   offset = sizestruct.offset;
   stride = sizestruct.stride;
   block  = sizestruct.block;
   count  = sizestruct.count; 
   ndim = length(offset);
   sid  = H5D.get_space(did);
   H5S.select_hyperslab(sid, 'H5S_SELECT_SET', ...
                        fliplr(offset) , ...
                        fliplr(stride) , ...
                        fliplr(count)  , ...
                        fliplr(block)  );

   mid = H5S.create_simple(ndim, fliplr(count), []);

   var = H5D.read(did,'H5ML_DEFAULT',mid, sid,'H5P_DEFAULT');
   H5S.close(sid);
   H5S.close(mid);
else
   sid = 'H5S_ALL';
   mid = 'H5S_ALL';
   var = H5D.read(did,'H5ML_DEFAULT',mid, sid,'H5P_DEFAULT');
end
H5D.close(did);

% remove blank spaces
if strcmp(class(var),'char') && blank
   var = deblank(var);
end

H5F.close(h5fid);

return
end
