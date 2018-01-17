
#include <iostream>
#include <cstdio>
#include "base_types.hpp"
#include "results.hpp"
#include "Readers/tiff.hpp"
#include "mpi.hpp"
#include "ui_calls.hpp"

namespace CCPi {

  void write_as_tiff(const std::string basename, const voxel_data &voxels,
		     const int offset, const unsigned int max_value,
		     const unsigned int width, const bool clamp);
  void write_real(const std::string basename, const voxel_data &voxels,
		  const int offset);
  void write_bgs(const std::string basename, const voxel_data &voxels,
		 const real voxel_origin[3], const real voxel_size[3],
		 const int offset, const int nz_voxels);

}

void CCPi::write_results(const std::string basename, const voxel_data &voxels,
			 const real voxel_origin[3], const real voxel_size[3],
			 const int offset, const int nz_voxels,
			 const output_format format, const bool clamp)
{
  switch (format) {
  case no_output:
    break;
  case unsigned_byte_tiff:
    write_as_tiff(basename, voxels, offset, 255, 8, clamp);
    break;
  case unsigned_short_tiff:
    write_as_tiff(basename, voxels, offset, 65535, 16, clamp);
    break;
  case native_dump:
    if (machine::get_number_of_processors() > 1)
      report_error("Format does not support distributed memory");
    else
      write_real(basename, voxels, offset);
    break;
  case signed_short_tiff:
    write_as_tiff(basename, voxels, offset, 32767, 16, clamp);
    break;
  case bgs_float_dump:
    if (machine::get_number_of_processors() > 1)
      report_error("Format does not support distributed memory");
    else
      write_bgs(basename, voxels, voxel_origin, voxel_size, offset, nz_voxels);
    break;
  default:
    report_error("Unknown output format");
    break;
  }
}

void CCPi::write_as_tiff(const std::string basename, const voxel_data &voxels,
			 const int offset, const unsigned int max_value,
			 const unsigned int width, const bool clamp)
{
  // Problem normalising output over range of data if we only have a
  // sub-block of the voxels.
  if (offset > 0)
    report_error("Tiff data range issue with blocks - Todo");
  if (width != 8 and width != 16)
    report_error("Width not supported for tiff writing");
  else {
    const voxel_data::size_type *s = voxels.shape();
	//Set the progress
	initialise_progress(s[2], "Saving data...");
    std::size_t n = s[0] * s[1];
    // copy buffer for tiff image slice
    unsigned short *sdata = new unsigned short[n];
    unsigned char *cdata = (unsigned char *)sdata;
    voxel_type vmax = -1e10;
    voxel_type vmin = +1e10;
    if (!clamp) {
      // find range to scale
      for (int i = 0; i < (int)s[0]; i++) {
	for (int j = 0; j < (int)s[1]; j++) {
	  for (int k = 0; k < (int)s[2]; k++) {
	    if (vmax < voxels[i][j][k])
	      vmax = voxels[i][j][k];
	    if (vmin > voxels[i][j][k])
	      vmin = voxels[i][j][k];
	  }
	}
      }
    }
    voxel_type scale = (voxel_type)max_value / (vmax - vmin);
    char index[8];
    bool ok = true;
    for (int k = 0; (k < (int)s[2] and ok); k++) {
      snprintf(index, 8, "_%04d", offset + k + 1);
      std::string name = basename + index + ".tif";
      int idx = 0;
      if (width == 8) {
	for (int j = 0; j < (int)s[1]; j++) {
	  for (int i = 0; i < (int)s[0]; i++) {
	    if (clamp) {
	      if (voxels[i][j][k] < 0.0)
		cdata[idx] = 0;
	      else if (voxels[i][j][k] >= 1.0)
		cdata[idx] = (unsigned char)max_value;
	      else
		cdata[idx] = (unsigned char) (voxels[i][j][k]
					      * (voxel_type)max_value);
	    } else
	      cdata[idx] = (unsigned char) ((voxels[i][j][k] - vmin) * scale);
	    idx++;
	  }
	}
      } else {
	for (int j = 0; j < (int)s[1]; j++) {
	  for (int i = 0; i < (int)s[0]; i++) {
	    if (clamp) {
	      if (voxels[i][j][k] < 0.0)
		sdata[idx] = 0;
	      else if (voxels[i][j][k] >= 1.0)
		sdata[idx] = (unsigned short)max_value;
	      else
		sdata[idx] = (unsigned short) (voxels[i][j][k]
					       * (voxel_type)max_value);
	    } else
	      sdata[idx] = (unsigned short) ((voxels[i][j][k] - vmin) * scale);
	    idx++;
	  }
	}
      }
      ok = write_tiff(name, cdata, (int)s[0], (int)s[1], width);
	  update_progress(k);
    }
	update_progress(s[2]);
    delete [] sdata;
  }
}

void CCPi::write_real(const std::string basename, const voxel_data &voxels,
		      const int offset)
{
  //std::cout << "start dump\n";
  std::string name = basename + ".dat";
  const voxel_data::size_type *s = voxels.shape();
  std::size_t n = s[0] * s[1] * s[2];
  std::FILE *file = fopen(name.c_str(), "a+b");
  if (file == 0)
    report_error(" Failed to open output file - ", name);
  else {
    std::size_t o = s[0] * s[1] * std::size_t(offset) * sizeof(voxel_type);
#ifdef WIN32
    _fseeki64(file, o, SEEK_SET);
#else
    fseek(file, o, SEEK_SET);
#endif // WIN32
    fwrite(voxels.data(), sizeof(voxel_type), n, file);
    fclose(file);
  }
  //std::cout << "end dump\n";
}

void CCPi::write_bgs(const std::string basename, const voxel_data &voxels,
		     const real voxel_origin[3], const real voxel_size[3],
		     const int offset, const int nz_voxels)
{
  // Assumes linear update of offset so each new block goes at end - Todo?
  //std::cout << "start dump\n";
  std::string name = basename + ".dat";
  const voxel_data::size_type *s = voxels.shape();
  std::size_t n = s[0] * s[1] * s[2];
  //Set the progress
  initialise_progress(s[2], "Saving data...");
  char mode[4];
  mode[3] = '\0';
  if (offset == 0) {
    mode[0] = 'w';
    mode[1] = 'b';
    mode[2] = '\0';
  } else {
    mode[0] = 'a';
    mode[1] = '+';
    mode[2] = 'b';
  }
  std::FILE *file = fopen(name.c_str(), mode);
  if (file == 0)
    report_error(" Failed to open output file - ", name);
  else {
    if (offset == 0) {
      fprintf(file, "%d %d %d\n", int(s[0]), int(s[1]), nz_voxels);
      // centre value in voxel rather than on edge, since not writing n+1
      // values to get final voxel boundary.
      real shift[3];
      shift[0] = voxel_origin[0] + voxel_size[0] / 2.0;
      shift[1] = voxel_origin[1] + voxel_size[1] / 2.0;
      shift[2] = voxel_origin[2] + voxel_size[2] / 2.0;
      // for k.p need to scale by bohr to angstrom
      fprintf(file, "%12.6f %12.6f %12.6f\n", shift[0], shift[1], shift[2]);
      fprintf(file, "%12.6f 0.0 0.0\n", voxel_size[0]);
      fprintf(file, "0.0 %12.6f 0.0\n", voxel_size[1]);
      fprintf(file, "0.0 0.0 %12.6f\n", voxel_size[2]);
      fprintf(file, "Image\n");
    }
    float *x = new float[n];
    sl_int l = 0;
    for (std::size_t k = 0; k < s[2]; k++) {
      for (std::size_t j = 0; j < s[1]; j++) {
	for (std::size_t i = 0; i < s[0]; i++) {
	  x[l] = (float)voxels[i][j][k];
	  l++;
	}
      }
	  update_progress(k);
	  fwrite(x+k*s[1]*s[0], sizeof(float), s[1] * s[0], file);
    }
	update_progress(s[2]);
    // would need use to store an offset for end of header block to work
    //std::size_t o = s[0] * s[1] * std::size_t(offset);
    //fseek(file, o, SEEK_CUR);
    //fseek(file, 0, SEEK_END); achieved by fopen("a")
    //fwrite(x, sizeof(float), n, file);
    delete [] x;
    fclose(file);
  }
  //std::cout << "end dump\n";
}
