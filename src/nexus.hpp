
#ifndef CCPI_RECON_NEXUS
#define CCPI_RECON_NEXUS

namespace CCPi {

  bool read_NeXus(pixel_data &pixels, pixel_2d &i_dark, pixel_2d &f_dark,
		  pixel_2d &i_bright, pixel_2d &f_bright, int &nh_pixels,
		  int &nv_pixels, std::vector<real> &angles, int &nangles,
		  real &hsize, real &vsize, const std::string filename,
		  const bool all_angles, const bool read_data,
		  const int start_idx, const int block_size);
  // Todo - write NeXus NXtomoproc?

}

#endif // CCPI_RECON_NEXUS
