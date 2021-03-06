function [] = write_float(w, pathname, filename, nvoxels)
% function to write reconstruction as raw float values

% input:
% x: voxel result
% pathname: name of path where files are stored
% filename: name of file to output

  % col-major [z y x] -> [x y z]
  data = permute(w,[3 2 1]);
  fid = fopen([pathname filename '.vol'], 'w');
  fwrite(fid, w, 'single');
  fclose(fid);

end
