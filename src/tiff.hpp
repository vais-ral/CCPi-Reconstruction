
#ifndef CCPI_RECON_TIFF
#define CCPI_RECON_TIFF

namespace CCPi {

  bool read_tiff(const std::string filename, pixel_type pixel_data[],
		 const int n_h_pixels, const int n_v_pixels);
  bool write_tiff(const std::string filename, unsigned short sdata[],
		  unsigned char cdata[], const int nx, const int ny,
		  const int width);

}

#endif // CCPI_RECON_TIFF