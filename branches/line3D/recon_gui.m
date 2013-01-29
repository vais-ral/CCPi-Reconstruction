% Initial attempt at a user interface to cgls reconstruction
% derived from recon3D
% BGS 24/8/12

[filename pathname] = uigetfile({'*.xtekct'}, 'Select XTek data');

addpath c/
addpath tools/

% add trailing directory separator to pathname
pathname = strcat(pathname, '/');
% strip .xtekct from filename to get base
filebase = strrep(filename, '.xtekct', '');

[nvoxels iterations beam_harden ok] = input_recon_data();

if ok

  voxels = [nvoxels nvoxels nvoxels];

% Load data and geometrical parameters from file
% Use this if you are using a reconstruction.xtekct file from FDK reconstruction
%[data geom] = load_data(pathname, filename, pathname2, filename2);
% Use this otherwise
  [data geom] = load_data(pathname, filebase);

% Find centre of rotation
  geom = centre_geom(data, geom);
  b = data(:);

% Optional beam-hardening correction
  if beam_harden
    b = b.^2;
  end

% Call cgls main function
  [cglsOut rho eta cancel] = cgls_XTek_single(b, iterations, geom, voxels);

  [m index] = min(rho);
  
  if index < size(rho,2)
      h = warndlg('Iterations diverged', 'CGLS', modal);
      uiwait(h);
  end
  
  if not(cancel)
    [filetype basename] = save_recon();

    if filetype == 1
      write_tiff(cglsOut(:,end), pathname, basename, voxels, 16);
    elseif filetype == 2
      write_tiff(cglsOut(:,end), pathname, basename, voxels, 8);
    else
      write_float(cglsOut(:,end), pathname, basename, voxels);
    end

  end

end
