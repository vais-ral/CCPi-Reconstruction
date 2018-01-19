
#ifndef CCPI_VOXEL_SETUP
#define CCPI_VOXEL_SETUP
#include "CCPiDefines.h"

CCPI_EXPORT void calculate_block_sizes(int &nx_voxels, int &ny_voxels,
				  int &nz_voxels, int &maxz_voxels,
				  int &block_size, int &block_step,
				  const int num_processors,
				  const int blocking_factor,
				  const int pixels_per_voxel,
				  CCPi::instrument *instrument,
				  const bool recon_blocks);

CCPI_EXPORT voxel_data *reconstruct(CCPi::instrument *device,
			       CCPi::reconstruction_alg *algorithm,
			       const std::string data_file,
			       const std::string output_name,
			       real full_vox_origin[3], real voxel_size[3],
			       const real rotation_centre,
			       const int pixels_per_voxel,
			       const int blocking_factor,
			       const bool beam_harden,
			       const CCPi::output_format write_format,
			       const bool clamp_output, const bool phantom);

CCPI_EXPORT voxel_data *reconstruct(CCPi::instrument *device,
			       CCPi::reconstruction_alg *algorithm,
			       const numpy_3d &pixels,
			       const numpy_1d &angles,
			       const real rotation_centre,
			       const int pixels_per_voxel,
			       const int blocking_factor,
			       const bool beam_hardenconst,const bool is_pixel_in_log);

CCPI_EXPORT voxel_data *calculate_dimensions(CCPi::instrument *device,
			       CCPi::reconstruction_alg *algorithm,
			       const numpy_3d &pixels,
			       const numpy_1d &angles,
			       const real rotation_centre,
			       const int pixels_per_voxel,
			       const int blocking_factor);

CCPI_EXPORT voxel_data *reconstruct(CCPi::instrument *device,
			       CCPi::reconstruction_alg *algorithm,
			       const numpy_3d &pixels,
			       const numpy_1d &angles,
			       const numpy_1d &h_offsets,
			       const numpy_1d &v_offsets,
			       const int pixels_per_voxel,
			       const real source_x, const real detector_x,
			       const real pixel_h_size, const real pixel_v_size,
			       const real mask_radius, const bool beam_harden,
			       real full_vox_origin[3], real voxel_size[3],
			       const bool has_offsets, const bool is_pixel_in_log);

/*Iteration with input solution at previous step*/
CCPI_EXPORT voxel_data *reconstruct(CCPi::instrument *device,
	CCPi::reconstruction_alg *algorithm,
	const numpy_3d &pixels, const numpy_1d &angles,
	const real rotation_centre, const int pixels_per_voxel,
	const int blocking_factor, const bool beam_harden, const bool is_pixel_in_log,
	voxel_data * last_iteration_voxels);	
#endif // CCPI_VOXEL_SETUP
