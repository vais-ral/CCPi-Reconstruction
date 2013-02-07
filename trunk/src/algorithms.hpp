
#ifndef CCPI_RECON_ALGORITHMS
#define CCPI_RECON_ALGORITHMS

namespace CCPi {

  enum algorithms { alg_FDK, alg_CGLS, alg_TVreg };

  bool cgls_reconstruction(const instrument *device, voxel_data &voxels,
			   const real origin[3], const real voxel_size[3]);

}

#endif // CCPI_RECON_ALGORITHMS
