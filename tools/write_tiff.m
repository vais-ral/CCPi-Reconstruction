function [] = write_tiff(w, pathname, filename, nvoxels, nbits)
% function to write reconstruction as a set of tiff images.

% input:
% x: voxel result
% pathname: name of path where files are stored
% filename: name of files to output
% nbits = 8 or 16 bit output

  x = w;
  dxmax = max(w);
  dxmin = min(w);
  % shift to 0.0
  x(:) = x(:) - dxmin;
  % scale to max of uint16
  %x(:) = x(:) * (65535 / dxmax);
  % AVS issue
  x(:) = x(:) * (32767 / dxmax);
  x = reshape(x,nvoxels);
  % write projection data
  for i = 1:nvoxels(3)
    data = x(:,:,i);
    if nbits == 8
      data = uint8(data);
    else
      data = uint16(data);
    end
    imwrite(data, [pathname filename '_' dec2base(i,10,4) '.tif']);
  end

end
