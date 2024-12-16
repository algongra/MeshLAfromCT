%  A simple example performing CPD using the wavefront OBJ file format.

% Clear figures, etc.
clear all; close all; clc;

% Paths to fixed and moving examples.
f_fn = 'C:\Users\woodford\ucsd_cvi_repo\tools\wobj\examples-dv\suzanne-fixed.obj';
m_fn = 'C:\Users\woodford\ucsd_cvi_repo\tools\wobj\examples-dv\suzanne-moving.obj';
o_fn = 'registered.obj';

% Define the CPD options.
opt.method='nonrigid_lowrank'; % Use the lowrank matrix approximation.
opt.numeig=30;                 % leave only 30 larges (out of 8171) eigenvectors/values to approximate G
opt.eigfgt=1;                  % use FGT to find the largest eigenvectore/values
opt.beta=2;                    % the width of Gaussian kernel (smoothness)
opt.lambda=3;                  % regularization weight
opt.viz=1;                     % show every iteration
opt.outliers=0.1;              % use 0.7 noise weight
opt.fgt=2;                     % use FGT to compute matrix-vector products (2 means to switch to truncated version at the end, see cpd_register)
opt.normalize=1;               % normalize to unit variance and zero mean before registering (default)
opt.corresp=0;                 % compute correspondence vector at the end of registration (not being estimated by default)
opt.max_it=100;                % max number of iterations
opt.tol=1e-3;                  % tolerance

cpd_register_obj(f_fn, m_fn, o_fn, opt)
