% Example 4. 3D Rigid CPD point-set registration. Full options intialization.
% 3D face point-set, outliers, and missing points.
clear all; close all; clc;

load cpd_data3D_face.mat; Y=X;

% add a random rigid transformation
R=cpd_R(rand(1),rand(1),rand(1));
X=rand(1)*X*R'+1;

% add outliers and delete some points
X=[X(1:end-20,:); 0.1*randn(40,3)]-2;
Y=[Y(20:end,:); 0.1*randn(40,3)]-2;

% Set the options
opt.method='rigid'; % use rigid registration
opt.viz=1;          % show every iteration
opt.outliers=0.6;   % use 0.6 noise weight to add robustness 

opt.normalize=1;    % normalize to unit variance and zero mean before registering (default)
opt.scale=1;        % estimate global scaling too (default)
opt.rot=1;          % estimate strictly rotational matrix (default)
opt.corresp=0;      % do not compute the correspondence vector at the end of registration (default)

opt.max_it=100;     % max number of iterations
opt.tol=1e-8;       % tolerance


% registering Y to X
[Transform, Correspondence]=cpd_register(X,Y,opt);


figure,cpd_plot_iter(X, Y); title('Before');

% X(Correspondence,:) corresponds to Y
figure,cpd_plot_iter(X, Transform.Y);  title('After registering Y to X');
