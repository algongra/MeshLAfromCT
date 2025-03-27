function [mesh,pp] = Mesh_SmoothSpline(mesh,dt_new,varargin);
% Get splines to reconstruct geometry
% OPTIONAL ARGUMENTS
% ------------------
%  - seg_time:    time instants at which segmentations were obtained along
%                 the cardiac cycle in [s]. [double].
%                 Size: 1xnt, where nt is the number of segmentations used
%                 ACHTUNG: seg_time MUST NOT have time instants from two
%                          different cardiac cycles (e.g., 0 and T MUST NOT be
%                          included, only one of those two values MUST be
%                          included). Default: []
%  - ifmodified:  Uses modified parameters defined by Yvonne
%  - ifcubicspln: Uses cubic spline interpolation (equivalent to use spline
%                 Matlab function to obtain splines)
% defaults
seg_time = [];
ifmodified = false;
ifcubicspln = false;
assigndefaults(varargin{:});
if isempty(seg_time)
    dt = mesh.dt;
    seg_time = 0:dt:1-dt;
end

aa = mesh.v_mov;
nt = mesh.nt;
nelem = size(aa,1);
ndim = size(aa,2);
n_tsteps_rep = nt - 1;


t = 0:dt_new:1; % where we want to evaluate

% Repeat some tsteps of the geometry at the beginning at the end to avoid
% problem with the periodicity
time = [-1+seg_time(end-n_tsteps_rep+1:end) ...
           seg_time ...
         1+seg_time(1:n_tsteps_rep)];
pts = cat(3,aa(:,:,end-n_tsteps_rep+1:end), ...
            aa(:,:,1:nt),...
            aa(:,:,1:n_tsteps_rep));

% Start off with the trapezoidal error weights
dt_seg = diff(time);
weights = ([dt_seg 0]+[0 dt_seg])/2; % weight is dt/2 at beginning and end

% csaps chooses the smoothing parameter (since p = -1)
if ~ifmodified
   p = -1;
   if ifcubicspln
      p = 1;
   end
   w = weights;
end % see below
xx = []; % No sampling (we will use fnval later)

% allocate memory
%pp.coefs = zeros(nelem,ndim,nt+1,4);
pp.coefs = zeros(nelem,ndim,length(time)-1,4);
pts = double(pts);

disp('Generating splines...')
tic
for e = 1:nelem
    for dim = 1:ndim
        if ifmodified
            n = abs(squeeze(pts(e,dim,2:end) - pts(e,dim,1:end-1))'...
                ./(time(2:end) - time(1:end-1)));
            n = max(n/max(n), 1e-10); % avoid zeros
            p = [-1 n.^(-1)];
            w = max(weights .* squeeze(abs(pts(e,dim,:) - mean(pts(e,dim,:),'all')))', 1e-10); % avoid zeros
        end
        spln = csaps(time,pts(e,dim,:),p,xx);
        %pp.coefs(e,dim,:,:) = spln.coefs(1+n_tsteps_rep:1+n_tsteps_rep+nt,:);
        pp.coefs(e,dim,:,:) = spln.coefs;
    end
end
tgenspln = toc;
disp(sprintf('elapsed time: %g',tgenspln));

%pp.breaks = spln.breaks(1+n_tsteps_rep:1+n_tsteps_rep+nt+1);
pp.breaks = spln.breaks;
pp.dim = spln.dim;
pp.form = spln.form;
pp.order = spln.order;
pp.pieces = size(pp.coefs,3);

aa = pp;
mesh.v_mov_old = mesh.v_mov;
mesh.v_mov_new = zeros(nelem,ndim,length(t));
disp('Evaluating splines...')
tic
for e = 1:nelem
    for dim = 1:ndim
        aa.coefs = squeeze(pp.coefs(e,dim,:,:));
        vspln = fnval(aa,t);
        mesh.v_mov_new(e,dim,:) = squeeze(vspln);
    end
end
tevalspln = toc;
disp(sprintf('elapsed time: %g',tevalspln));
mesh.v_mov = mesh.v_mov_new(:,:,1:end-1);
mesh.dt = dt_new;
mesh.nt = size(mesh.v_mov,3);

fprintf('\n\nerror in spline periodicity = %g, %g, %g\n',...
    max(abs(mesh.v_mov_new(:,:,1)-mesh.v_mov_new(:,:,end)),[],1));
%duration(datetime('now','Format','mm:ss.SSS')-tic,'Format','mm:ss.SSS')
if isfield(mesh,'normals_mov')
    mesh = rmfield(mesh,{'normals_mov', 'areas_mov'});
end

% Get mesh properties (area,quality and normals)
mesh = Mesh_AnalyzePM(mesh);
