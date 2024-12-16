function [xyz da vol] = points_on_sphere_leopardi(D,dx)
%---------------------------------------------------------------------%
% determines an "even" distribution of points *on*the*surface* of a sphere
% the data is written to ascii file in a 4-column format (x,y,z)
%
% this code calls the function "eq_point_set" from the package "EQSP" by
% paul leopardi
% http://www.mathworks.com/matlabcentral/fileexchange/13356-eqsp-recursive-zonal-sphere-partitioning-toolbox
%---------------------------------------------------------------------%
path(path,'/data/mmoriche/IBcode/_unregistered/thirdparty/eq_sphere_partitions/eq_utilities');
path(path,'/data/mmoriche/IBcode/_unregistered/thirdparty/eq_sphere_partitions/eq_test');
path(path,'/data/mmoriche/IBcode/_unregistered/thirdparty/eq_sphere_partitions/eq_region_props');
path(path,'/data/mmoriche/IBcode/_unregistered/thirdparty/eq_sphere_partitions/eq_point_set_props');
path(path,'/data/mmoriche/IBcode/_unregistered/thirdparty/eq_sphere_partitions/eq_partitions');
path(path,'/data/mmoriche/IBcode/_unregistered/thirdparty/eq_sphere_partitions/eq_illustrations');
path(path,'/data/mmoriche/IBcode/_unregistered/thirdparty/eq_sphere_partitions');
%%---------------------------------------------------------------------%
%points_per_diameter=72;
%
radius=D/2;
%points_per_diameter=48;
points_per_diameter=round(2*radius/dx);
outer_radius=radius+dx/2;
points=[];
%
lrad=radius;
np=floor((pi/3)*(12*lrad^2/dx^2+1));
nls=np;
tmp=eq_point_set(2,np); %<-this is from the "EQ_POINT_SET" package
tmp=tmp*radius;    %<-scale to present radius
dvol=pi*dx/(3*np)*(12*lrad^2+dx^2); %shell between lrad(+/-)dx/2
%tmp(4,:)=ones(size(tmp,2),1)*dvol;
points=tmp;
disp(sprintf('total number of points: %d',nls))
%
%----------------------------------------------------------------
%%figure;plot3(points(1,:),points(2,:),points(3,:),'ko');axis image;
xyz = points';
da = 4*pi*radius^2/nls;
vol = dvol;
%----------------------------------------------------------------
%
% how to judge the quality of the distribution?
%
% 1. average point position???
%
return
end
