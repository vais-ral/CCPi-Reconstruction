function [b geom] = create_phantom(test)
% Function to create simple cube phantom for testing the CGLS
% reconstruction code

% W. Thompson (adapted by N. Wadeson)

% 03/04/2012

disp('Creating phantom...') 
% set up geometry (realistic values for small XTek machine)
geom.source.x = -250;
geom.source.y = 0;
geom.source.z = 0;

geom.dets.x = 737;


if (nargin > 0) % set fast test parameters
    geom.dets.y = (-116.5225:1:116.5225)';
    geom.dets.z = (-92.3925:1:92.3925)';
 
   % set up grid (unrealistic resolution for test)
    voxels = [10 10 10];
else % set realistic parameters
    geom.dets.y = (-116.5225:0.127:116.5225)';
    geom.dets.z = (-92.3925:0.127:92.3925)';

    % set up grid (use low resolution to save time)
    voxels = [459 459 359];
end

geom.d_sd = 987;

geom.angles = linspace(0,2*pi,501)';
geom.angles = geom.angles(1:500);

halfDetSize = 0.5*(geom.dets.y(end) - geom.dets.y(end-1));
geom.mask_radius = (-geom.source.x) * sin(atan((geom.dets.y(end)+halfDetSize)/geom.d_sd));
geom.voxel_size = (2*geom.mask_radius/voxels(1))*[1 1 1];
image_vol = voxels.*geom.voxel_size;
image_offset = -image_vol/2;

% set up phantom volume
x = single(zeros(voxels));

if (nargin > 0) % set fast test parameters
    % add cubes to phantom volume
    x(1:5,1:5,1:5) = 1;
    x(6:10,6:10,6:10) = 1;
else
    % add cubes to phantom volume
    x(108:189,108:189,58:139) = 1;
    x(190:271,190:271,140:221) = 1;
    x(272:353,272:353,222:303) = 1;
end

x = x(:);

tic
% perform projection step
b = CBproject_single(x, geom, voxels, geom.voxel_size, image_offset);
toc

% add noise
noise_level = 0.01;
noise = randn(length(b),1);
b = b + noise_level*std(b)*noise;

disp('Phantom created.') 


