
#ifndef CCPI_VOXEL_SETUP
#define CCPI_VOXEL_SETUP

extern void calculate_block_sizes(int &nx_voxels, int &ny_voxels,
				  int &nz_voxels, int &maxz_voxels,
				  int &block_size, int &block_step,
				  const int num_processors,
				  const int blocking_factor,
				  const int pixels_per_voxel,
				  CCPi::instrument *instrument,
				  const bool recon_blocks);

extern void reconstruct(CCPi::instrument *device,
			CCPi::reconstruction_alg *algorithm,
			const std::string data_file,
			const std::string output_name,
			const real rotation_centre, const int pixels_per_voxel,
			const int blocking_factor, const bool beam_harden,
			const CCPi::output_format write_format,
			const bool clamp_output, const bool phantom);

#endif // CCPI_VOXEL_SETUP
