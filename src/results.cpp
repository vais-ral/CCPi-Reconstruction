
#include <iostream>
#include <cstdio>
#include "base_types.hpp"
#include "results.hpp"
#include "tiff.hpp"

namespace CCPi {

  void write_as_tiff(const std::string basename, const voxel_data &voxels,
		     const unsigned int max_value, const unsigned int width);
  void write_real(const std::string basename, const voxel_data &voxels);
  void write_bgs(const std::string basename, const voxel_data &voxels,
		 const real voxel_origin[3], const real voxel_size[3]);

}

void CCPi::write_results(const std::string basename, const voxel_data &voxels,
			 const real voxel_origin[3], const real voxel_size[3],
			 const output_format format)
{
  switch (format) {
  case unsigned_byte_tiff:
    write_as_tiff(basename, voxels, 255, 8);
    break;
  case unsigned_short_tiff:
    write_as_tiff(basename, voxels, 65535, 16);
    break;
  case native_dump:
    write_real(basename, voxels);
    break;
  case signed_short_tiff:
    write_as_tiff(basename, voxels, 32767, 16);
    break;
  case bgs_float_dump:
    write_bgs(basename, voxels, voxel_origin, voxel_size);
    break;
  default:
    std::cerr << "Unknown output format\n";
    break;
  }
}

void CCPi::write_as_tiff(const std::string basename, const voxel_data &voxels,
			 const unsigned int max_value, const unsigned int width)
{
  if (width != 8 and width != 16)
    std::cerr << "Width not supported for tiff writing\n";
  else {
    const voxel_data::size_type *s = voxels.shape();
    std::size_t n = s[0] * s[1];
    // copy buffer for tiff image slice
    unsigned short *sdata = new unsigned short[n];
    unsigned char *cdata = (unsigned char *)sdata;
    // find range to scale - or should we truncate on 0.0/1.0?
    voxel_type vmax = -1e10;
    voxel_type vmin = +1e10;
    for (int k = 0; k < (int)s[2]; k++) {
      for (int j = 0; j < (int)s[1]; j++) {
	for (int i = 0; i < (int)s[0]; i++) {
	  if (vmax < voxels[i][j][k])
	    vmax = voxels[i][j][k];
	  if (vmin > voxels[i][j][k])
	    vmin = voxels[i][j][k];
	}
      }
    }
    voxel_type scale = (voxel_type)max_value / (vmax - vmin);
    char index[8];
    bool ok = true;
    for (int k = 0; (k < (int)s[2] and ok); k++) {
      snprintf(index, 8, "_%04d", k + 1);
      std::string name = basename + index + ".tif";
      int idx = 0;
      if (width == 8) {
	for (int j = 0; j < (int)s[1]; j++) {
	  for (int i = 0; i < (int)s[0]; i++) {
	    /*
	    if (voxels[i][j][k] < 0.0)
	      cdata[idx] = 0;
	    else if (voxels[i][j][k] > 1.0)
	      cdata[idx] = (unsigned char)max_value;
	    else
	      cdata[idx] = (unsigned char) (voxels[i][j][k]
					     * (voxel_type)max_value);
	    */
	    cdata[idx] = (unsigned char) ((voxels[i][j][k] - vmin) * scale);
	    idx++;
	  }
	}
      } else {
	for (int j = 0; j < (int)s[1]; j++) {
	  for (int i = 0; i < (int)s[0]; i++) {
	    /* This would truncate
	    if (voxels[i][j][k] < 0.0)
	      sdata[idx] = 0;
	    else if (voxels[i][j][k] > 1.0)
	      sdata[idx] = (unsigned short)max_value;
	    else
	      sdata[idx] = (unsigned short) (voxels[i][j][k]
					     * (voxel_type)max_value);
	    */
	    sdata[idx] = (unsigned short) ((voxels[i][j][k] - vmin) * scale);
	    idx++;
	  }
	}
      }
      ok = write_tiff(name, sdata, cdata, (int)s[0], (int)s[1], width);
    }
    delete [] sdata;
  }
}

void CCPi::write_real(const std::string basename, const voxel_data &voxels)
{
  std::cout << "start dump\n";
  std::string name = basename + ".dat";
  const voxel_data::size_type *s = voxels.shape();
  std::size_t n = s[0] * s[1] * s[2];
  std::FILE *file = fopen(name.c_str(), "w");
  if (file == 0)
    std::cerr << " Failed to open output file - " << name << '\n';
  else {
    fwrite(voxels.data(), sizeof(voxel_type), n, file);
    fclose(file);
  }
  std::cout << "end dump\n";
}

void CCPi::write_bgs(const std::string basename, const voxel_data &voxels,
		     const real voxel_origin[3], const real voxel_size[3])
{
  std::cout << "start dump\n";
  std::string name = basename + ".dat";
  const voxel_data::size_type *s = voxels.shape();
  std::size_t n = s[0] * s[1] * s[2];
  std::FILE *file = fopen(name.c_str(), "w");
  if (file == 0)
    std::cerr << " Failed to open output file - " << name << '\n';
  else {
    float *x = new float[n];
    for (std::size_t i = 0; i < n; i++)
      x[i] = (float)((voxels.data())[i]);
    fprintf(file, "%d %d %d\n", int(s[0]), int(s[1]), int(s[2]));
    // centre value in voxel rather than on edge, since not writing n+1
    // values to get final voxel boundary.
    real shift[3];
    shift[0] = voxel_origin[0] + voxel_size[0] / 2.0;
    shift[1] = voxel_origin[1] + voxel_size[1] / 2.0;
    shift[2] = voxel_origin[2] + voxel_size[2] / 2.0;
    // for k.p need to scale by bohr to angstrom
    fprintf(file, "%12.6f %12.6f %12.6f\n", shift[0] / 0.52917721092,
	    shift[1] / 0.52917721092, shift[2] / 0.52917721092);
    fprintf(file, "%12.6f 0.0 0.0\n", voxel_size[0] / 0.52917721092);
    fprintf(file, "0.0 %12.6f 0.0\n", voxel_size[1] / 0.52917721092);
    fprintf(file, "0.0 0.0 %12.6f\n", voxel_size[2] / 0.52917721092);
    fprintf(file, "Image\n");
    fwrite(x, sizeof(float), n, file);
    delete [] x;
    fclose(file);
  }
  std::cout << "end dump\n";
}
