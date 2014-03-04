
#ifndef CCPI_RECON_ALGORITHMS
#define CCPI_RECON_ALGORITHMS

#include <list>

namespace CCPi {

  enum algorithms { alg_FDK, alg_CGLS, alg_TVreg };

  bool fbp_reconstruction(const instrument *device, voxel_data &voxels,
			  const real origin[3], const real voxel_size[3],
			  const filter_name_t name,
			  const filter_window_t window,
			  const filter_norm_t norm, const real bandwidth);
  bool cgls_reconstruction(const class instrument *device, voxel_data &voxels,
			   const real origin[3], const real voxel_size[3],
			   const int iterations);
  bool tv_regularization(const instrument *device, voxel_data &voxels,
			 const real origin[3], const real voxel_size[3],
			 const real alpha, const real tau, const real init_L,
			 const real init_mu, const int constraint);

  void tvreg_core(voxel_type *xkp1, real *fxkp1, real *hxkp1, real *gxkp1,
		  real *fxkp1l, int *kend, const real voxel_size[],
		  const real *b, const real alpha, real tau,
		  real bL, real bmu, real epsb_rel,int k_max, const int Ddim,
		  const int Dm, const int Dn, const int Dl,
		  const sl_int prodDims, int ctype, real *d, real *c,
		  const bool ghxl, const bool xl, real *hxkp1l,
		  real *gxkp1l, real *xlist, const bool verbose,
		  real *numGrad, real* numBack, real *numFunc,
		  real *numRest, real *Lklist, real *muklist,
		  std::list<int> &rp, const real grid_offset[],
		  const class instrument *device);

}

#endif // CCPI_RECON_ALGORITHMS
