
#ifndef CCPI_RECON_ALGORITHMS
#define CCPI_RECON_ALGORITHMS

#include <list>

namespace CCPi {

  enum algorithms { alg_FDK, alg_CGLS, alg_TVreg, alg_BiCGLS, alg_BiCGSTABLS,
		    alg_landweber, alg_MLEM, alg_SIRT, alg_CGLS_Tikhonov };

  class reconstruction_alg {
  public:
    virtual ~reconstruction_alg();
    virtual bool reconstruct(class instrument *device, voxel_data &voxels,
			     const real origin[3],
			     const real voxel_size[3]) = 0;
    virtual bool supports_blocks() const = 0;
  };

}

inline CCPi::reconstruction_alg::~reconstruction_alg()
{
}

#endif // CCPI_RECON_ALGORITHMS
