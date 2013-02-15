
#ifndef CCPI_RECON_RESULTS
#define CCPI_RECON_RESULTS

namespace CCPi {

  enum output_format { signed_byte_tiff, unsigned_byte_tiff, signed_short_tiff,
		       unsigned_short_tiff, native_dump, bgs_float_dump };

  void write_results(const std::string basename, const voxel_data &voxels,
		     const real voxel_origin[3], const real voxel_size[3],
		     const output_format format);

}

#endif // CCPI_RECON_RESULTS
