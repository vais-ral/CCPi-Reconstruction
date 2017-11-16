function [x rho eta cancel] = cgls_XTek_single(b, iterations, geom, voxels)
% CGLS reconstruction code for parallel beam geometries

disp('Reconstructing...')

% check input b is single precision
b = single(b);

% set up reconstruction grid with default values for small XTek if not
% given as input
if nargin < 4
    voxels = [1836 1836 1436];
end

nx = (geom.dets.y(numel(geom.dets.y)) - geom.dets.y(1)) / (voxels(1) - 1);
ny = (geom.dets.y(numel(geom.dets.y)) - geom.dets.y(1)) / (voxels(2) - 1);
nz = (geom.dets.z(numel(geom.dets.z)) - geom.dets.z(1)) / (voxels(3) - 1);
geom.voxel_size = [nx ny nz];
image_offset = [geom.dets.y(1) geom.dets.y(1) geom.dets.z(1)] - geom.voxel_size/2;

n_vox = prod(voxels);

x = single(zeros(n_vox,1)); % storage for solution

rho = zeros(iterations,1);
eta = zeros(iterations,1);

outputIts = 10; % output solution after every 10 iterations

% with storage
x_stored = repmat(x,1,floor(iterations/outputIts));

% x_stored = repmat(x,1,iterations);

% Prepare for CG iteration.
h = waitbar(0,'Preparing for CG iteration...', 'CreateCancelBtn',...
           'setappdata(gcbf,''canceling'',1)' );
setappdata(h,'canceling',0);
disp('Preparing for CG iteration...')
cg_stopped = 0;

tic
d = PBbackproject_single(b, geom, voxels, geom.voxel_size, image_offset);
% d(vox_ind) = 0;
r = b; 
normr2 = d'*d;
toc

if getappdata(h,'canceling')
  cg_stopped = 1;
end

if cg_stopped == 0

  waitbar(0.5/(iterations+0.5),h,'CG Iterating...');

  % Iterate.
  for i = 1:iterations
     if getappdata(h,'canceling')
        cg_stopped = 1;
        break
      end
      disp(i)
      tic
      % Update x and r vectors. 
      Ad = PBproject_single(d, geom, voxels, geom.voxel_size, image_offset);
      alpha = normr2/(Ad'*Ad); 
      x  = x + alpha*d; 
      r  = r - alpha*Ad; 
      s  = PBbackproject_single(r, geom, voxels, geom.voxel_size, image_offset);

      % Update d vector. 
      normr2_new = s'*s; 
      beta = normr2_new/normr2; 
      normr2 = normr2_new; 
      d = s + beta*d;
      toc

      % Storing data for output
      if (floor(i/outputIts) == i/outputIts)
          x_stored(:,i/outputIts) = x;
      elseif (i == iterations)
          x_stored(:,end+1) = x;
      end

      waitbar((i+0.5)/(iterations+0.5), h);

      if (nargout>1), rho(i) = norm(r); end
      if (nargout>2), eta(i) = norm(x); end
  end

end

delete(h)

x = x_stored;

if cg_stopped == 0
  disp('Reconstruction successfully completed.')
  cancel = false;
else
  disp('Reconstruction aborted')
  cancel = true;
end
