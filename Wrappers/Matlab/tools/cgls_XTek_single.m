function [x rho eta cancel] = cgls_XTek_single(b, iterations, geom, voxels)
% CGLS reconstruction code for the XTek circular scan geometry

% 22/08/2011

% 08/09/2011 single precision version using C routines for the cone beam
% projection and back projection

% 11/01/2012 added optional input argument for dimensions of the
% reconstruction volume in voxels

disp('Reconstructing...')

% check input b is single precision
b = single(b);

% set up reconstruction grid with default values for small XTek if not
% given as input
if nargin < 4
    voxels = [1836 1836 1436];
end

if (length(fieldnames(geom)) == 5) % data from file (only one file loaded)
    %maxDets_y = 2000;
    %geom.voxel_size = ((2*geom.mask_radius)/maxDets_y)*[1 1 1];
    %voxel_size = (2*(geom.mask_radius*(geom.dets.ny/maxDets_y))/voxels(1))*[1 1 1];
    geom.voxel_size = (2*geom.mask_radius/voxels(1))*[1 1 1];
    geom.offset = [0.0 0.0 0.0];
elseif (length(fieldnames(geom)) == 6) % using phantom data
    geom.offset = [0.0 0.0 0.0];
elseif (length(fieldnames(geom)) == 8) % data from file (plus 2nd FDK reconstruction file)
    if (voxels(3) == 1)
        geom.nVoxels(3) = 1;
        geom.offset(3) = 0.0;


%  % ****COMMENT/UNCOMMENT*************
%  % For TV_reg cut-down only. Uncomment if you wish to reconstruct using
%  % FDK reconstruction file but with a smaller number of z voxels.
%  % Do the same in tvreg_upn.m
%      else
%          geom.nVoxels(3) = 41;
%  % ****COMMENT/UNCOMMENT*************
      end


    voxels = geom.nVoxels;
else
    disp('geom structures contains %d fields', length(fieldnames(geom)));
    return
end

image_vol = voxels.*geom.voxel_size;
image_offset = -image_vol/2 + geom.offset;

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
d = CBbackproject_single(b, geom, voxels, geom.voxel_size, image_offset);
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
      Ad = CBproject_single(d, geom, voxels, geom.voxel_size, image_offset);
      alpha = normr2/(Ad'*Ad); 
      x  = x + alpha*d; 
      r  = r - alpha*Ad; 
      s  = CBbackproject_single(r, geom, voxels, geom.voxel_size, image_offset);

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
