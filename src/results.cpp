
#include <iostream>
#include <cstdio>
#include "base_types.hpp"
#include "results.hpp"

namespace CCPi {

  void write_real(const std::string basename, const voxel_data &voxels);
  void write_bgs(const std::string basename, const voxel_data &voxels,
		 const real voxel_origin[3], const real voxel_size[3]);

}

void CCPi::write_results(const std::string basename, const voxel_data &voxels,
			 const real voxel_origin[3], const real voxel_size[3],
			 const output_format format)
{
  switch (format) {
  case signed_byte_tiff:
    std::cerr << "Todo\n";
    break;
  case unsigned_byte_tiff:
    std::cerr << "Todo\n";
    break;
  case signed_short_tiff:
    std::cerr << "Todo\n";
    break;
  case unsigned_short_tiff:
    std::cerr << "Todo\n";
    break;
  case native_dump:
    write_real(basename, voxels);
    break;
  case bgs_float_dump:
    write_bgs(basename, voxels, voxel_origin, voxel_size);
    break;
  default:
    std::cerr << "Unknown output format\n";
    break;
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
    // for k.p need to scale by bohr to angstrom
    fprintf(file, "%12.6f %12.6f %12.6f\n", voxel_origin[0] / 0.52917721092,
	    voxel_origin[1] / 0.52917721092, voxel_origin[2] / 0.52917721092);
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
