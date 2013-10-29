function [] = write_float(w, pathname, filename, nvoxels, voxel_size)
% function to write reconstruction as raw float values

% input:
% x: voxel result
% pathname: name of path where files are stored
% filename: name of file to output

  fid = fopen([pathname filename '.vol'], 'w');
  image_vol = voxel_size.*nvoxels;
  image_offset = - image_vol / 2.0;

  shift = image_offset + voxel_size / 2.0;
  fprintf(fid, '%d %d %d\n', nvoxels(1), nvoxels(2), nvoxels(3));
  fprintf(fid, '%12.6f %12.6f %12.6f\n', shift(1), shift(2), shift(3));
  fprintf(fid, '%12.6f 0.0 0.0\n', voxel_size(1));
  fprintf(fid, '0.0 %12.6f 0.0\n', voxel_size(2));
  fprintf(fid, '0.0 0.0 %12.6f\n', voxel_size(3));
  fprintf(fid, 'Image\n');
  fwrite(fid, w, 'single');
  fclose(fid);

end
