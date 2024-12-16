function data = CT_nii2vol(seg_fullfile, isotropic, fast_opt)

if nargin < 2
    isotropic = false;
end
if nargin < 3
    fast_opt = false;
end

nii_data = load_nii(seg_fullfile);

%% Get dimension and build rectangular grid
data.nx = nii_data.hdr.dime.dim(2);
data.ny = nii_data.hdr.dime.dim(3);
data.nz = nii_data.hdr.dime.dim(4);

data.dx = nii_data.hdr.dime.pixdim(2);
data.dy = nii_data.hdr.dime.pixdim(3);
data.dz = nii_data.hdr.dime.pixdim(4);

data.x = 0.5*data.dx:data.dx:(data.nx-0.5)*data.dx;
data.y = 0.5*data.dy:data.dy:(data.ny-0.5)*data.dy;
data.z = 0.5*data.dz:data.dz:(data.nz-0.5)*data.dz;

if isotropic
    ori_grid = data;
    data.dx = min([data.dx, data.dy, data.dz]);
    data.dy = data.dx;
    data.dz = data.dx;

    data.x = single(data.x(1):data.dx:data.x(end));
    data.y = single(data.y(1):data.dy:data.y(end));
    data.z = single(data.z(1):data.dz:data.z(end));

    data.nx = length(data.x);
    data.ny = length(data.y);
    data.nz = length(data.z);

    % The X (first dimension) coordinate is reversed
    img = single(flip(nii_data.img,1));

    if data.nz ~= ori_grid.nz
       data.img = permute(int8(interp1(ori_grid.z, permute(img,[3 2 1]), data.z, 'nearest')),[3 2 1]);
    else
       data.img = img;
    end
    if data.nx ~= ori_grid.nx || data.ny ~= ori_grid.ny
        error('Interpolation to isotropic grid implemented only across z');
    end
else
    % The X (first dimension) coordinate is reversed
    data.img = flip(nii_data.img,1);
end

if ~fast_opt
    % Y is the second dimension
    [data.X, data.Y,data.Z] = meshgrid(data.x,data.y,data.z);
end

% AGG Check (data.img and data.[X|Y|Z] are permute([2 1 3])
if sum(size(data.X) == size(data.img)) < 3
   data.img = permute(data.img,[2 1 3]);
end
