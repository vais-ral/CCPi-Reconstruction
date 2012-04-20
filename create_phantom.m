function [b geom] = create_phantom
% Function to create simple cube phantom for testing the CGLS
% reconstruction code

% W. Thompson

% 03/04/2012

% set up geometry (realistic values for small XTek machine)
geom.source.x = -250;
geom.source.y = 0;
geom.source.z = 0;

geom.dets.x = 737;
geom.dets.y = (-116.5225:0.127:116.5225)';
geom.dets.z = (-92.3925:0.127:92.3925)';

geom.d_sd = 987;

geom.angles = linspace(0,2*pi,501)';
geom.angles = geom.angles(1:500);

% set up grid (use low resolution to save time)
voxels = [459 459 359];

mask_radius = (-geom.source.x) * sin(atan(geom.dets.y(end)/geom.d_sd));
voxel_size = (2*mask_radius/voxels(1))*[1 1 1];
image_vol = voxels.*voxel_size;
image_offset = -image_vol/2;

% set up phantom volume
x = single(zeros(voxels));

% add cubes
x(108:189,108:189,58:139) = 1;
x(190:271,190:271,140:221) = 1;
x(272:353,272:353,222:303) = 1;

x = x(:);

% perform projection step
b = CBproject_single(x, geom, voxels, voxel_size, image_offset);

% add noise
noise_level = 0.01;
noise = randn(length(b),1);
b = b + noise_level*std(b)*noise;



