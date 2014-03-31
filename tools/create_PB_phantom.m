function [b geom] = create_PB_phantom(test)
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
else % set test parameters
    geom.dets.y = (-116.5225:0.508:116.5225)';
    geom.dets.z = (-92.3925:0.508:92.3925)';
    %geom.dets.z = (-0.254:0.508:0.254)';

    % set up grid (use low resolution to save time)
    voxels = [459 459 364];
    %voxels = [459 459 2];
end

geom.angles = linspace(0,pi,251)';
geom.angles = geom.angles(1:250);

geom.voxel_size = [(geom.dets.y(2) - geom.dets.y(1)) (geom.dets.y(2) - geom.dets.y(1)) (geom.dets.z(2) - geom.dets.z(1))];
image_vol = voxels.*geom.voxel_size;
image_offset = [geom.dets.y(1) geom.dets.y(1) geom.dets.z(1)] - geom.voxel_size/2;

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
    %x(108:189,108:189,1:1) = 1;
end

x = x(:);

tic
% perform projection step
b = PBproject_single(x, geom, voxels, geom.voxel_size, image_offset);
toc

% add noise
noise_level = 0.01;
noise = randn(length(b),1);
b = b + noise_level*std(b)*noise;

disp('Phantom created.') 


