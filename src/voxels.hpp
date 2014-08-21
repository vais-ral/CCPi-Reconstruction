
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

#endif // CCPI_VOXEL_SETUP
