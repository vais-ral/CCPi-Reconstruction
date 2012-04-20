function [x rho eta] = cgls_XTek_single(b, iterations, geom, voxels)
% CGLS reconstruction code for the XTek circular scan geometry

% 22/08/2011

% 08/09/2011 single precision version using C routines for the cone beam
% projection and back projection

% 11/01/2012 added optional input argument for dimensions of the
% reconstruction volume in voxels

% check input b is single precision
b = single(b);

% set up reconstruction grid with default values for small XTek if not
% given as input
if nargin < 4
    voxels = [1836 1836 1436];
end

mask_radius = (-geom.source.x) * sin(atan(geom.dets.y(end)/geom.d_sd));
voxel_size = (2*mask_radius/voxels(1))*[1 1 1];
image_vol = voxels.*voxel_size;
image_offset = -image_vol/2;

n_vox = prod(voxels);

x = single(zeros(n_vox,1)); % storage for solution

rho = zeros(iterations,1);
eta = zeros(iterations,1);

% with storage
x_stored = repmat(x,1,floor(iterations/10));
% x_stored = repmat(x,1,iterations);

% Prepare for CG iteration.
disp('Preparing for CG iteration...')
tic
d = CBbackproject_single(b, geom, voxels, voxel_size, image_offset);
% d(vox_ind) = 0;
r = b; 
normr2 = d'*d;
toc

% Iterate.
for i = 1:iterations
    disp(i)
    tic
    % Update x and r vectors. 
    Ad = CBproject_single(d, geom, voxels, voxel_size, image_offset);
    alpha = normr2/(Ad'*Ad); 
    x  = x + alpha*d; 
    r  = r - alpha*Ad; 
    s  = CBbackproject_single(r, geom, voxels, voxel_size, image_offset);


    % Update d vector. 
    normr2_new = s'*s; 
    beta = normr2_new/normr2; 
    normr2 = normr2_new; 
    d = s + beta*d;
    toc
        
    if floor(i/10) == i/10
        x_stored(:,i/10) = x;
    end
%     x_stored(:,i) = x;
    
    if (nargout>1), rho(i) = norm(r); end
    if (nargout>2), eta(i) = norm(x); end
end
x = x_stored;













