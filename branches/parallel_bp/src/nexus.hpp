
#ifndef CCPI_RECON_NEXUS
#define CCPI_RECON_NEXUS

namespace CCPi {

  bool read_NeXus(pixel_type * &pixels, pixel_type * &i_dark,
		  pixel_type * &f_dark, pixel_type * &i_bright,
		  pixel_type * &f_bright, int &nh_pixels, int &nv_pixels,
		  real * &angles, int &nangles, real &hsize, real &vsize,
		  const std::string filename, const bool all_angles);
  // Todo - write NeXus NXtomoproc?

}

#endif // CCPI_RECON_NEXUS
