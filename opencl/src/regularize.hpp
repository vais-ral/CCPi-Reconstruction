
#ifndef CCPI_REGULARIZERS
#define CCPI_REGULARIZERS

namespace CCPi {

  void tikhonov_regularize(voxel_data &b, const voxel_data &a,
			   const int nx, const int ny, const int nz);
  void tv_regularize(voxel_data &b, const voxel_data &a,
		     const int nx, const int ny, const int nz);

}

#endif // CCPI_REGULARIZERS
