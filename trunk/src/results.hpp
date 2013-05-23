
#ifndef CCPI_RECON_RESULTS
#define CCPI_RECON_RESULTS

namespace CCPi {

  enum output_format { unsigned_byte_tiff, unsigned_short_tiff,
		       native_dump, signed_short_tiff, bgs_float_dump };

  void write_results(const std::string basename, const voxel_data &voxels,
		     const real voxel_origin[3], const real voxel_size[3],
		     const output_format format, const bool clamp);

}

#endif // CCPI_RECON_RESULTS
