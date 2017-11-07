
#ifndef CCPI_RECON_TIFF
#define CCPI_RECON_TIFF

namespace CCPi {

  bool read_tiff(const std::string filename, pixel_data &pixels,
		 const int angle, const int n_h_pixels, const int n_v_pixels,
		 const int v_shift);
  bool write_tiff(const std::string filename, unsigned char data[],
		  const int nx, const int ny, const int width);

}

#endif // CCPI_RECON_TIFF
