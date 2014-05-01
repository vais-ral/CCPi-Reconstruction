
#ifndef CCPI_RECON_ALGORITHMS
#define CCPI_RECON_ALGORITHMS

#include <list>

namespace CCPi {

  enum algorithms { alg_FDK, alg_CGLS, alg_TVreg };

  class reconstruction_alg {
  public:
    virtual bool reconstruct(const class instrument *device, voxel_data &voxels,
			     const real origin[3],
			     const real voxel_size[3]) = 0;
  };

  class cgls_base : public reconstruction_alg {
  public:
    cgls_base(const int niterations);

    bool reconstruct(const class instrument *device, voxel_data &voxels,
		     const real origin[3], const real voxel_size[3]);

  private:
    int iterations;
  };

  class tv_regularization : public reconstruction_alg {
  public:
    tv_regularization(const real alph, const real t, const real L,
		      const real mu, const int c);

    bool reconstruct(const instrument *device, voxel_data &voxels,
		     const real origin[3], const real voxel_size[3]);

    static void tvreg_core(voxel_type *xkp1, real *fxkp1, real *hxkp1,
			   real *gxkp1, real *fxkp1l, int *kend,
			   const real voxel_size[], const pixel_type *b,
			   const real alpha, real tau, real bL, real bmu,
			   real epsb_rel,int k_max, const int Ddim,
			   const int Dm, const int Dn, const int Dl,
			   const sl_int prodDims, int ctype, real *d, real *c,
			   const bool ghxl, const bool xl, real *hxkp1l,
			   real *gxkp1l, real *xlist, const bool verbose,
			   real *numGrad, real* numBack, real *numFunc,
			   real *numRest, real *Lklist, real *muklist,
			   std::list<int> &rp, const real grid_offset[],
			   const class instrument *device);

  private:
    real alpha;
    real tau;
    real init_L;
    real init_mu;
    int constraint;

    bool reconstruct(const instrument *device, pixel_type *b,
		     voxel_data &voxels, const real origin[3],
		     const real voxel_size[3]);
  };

}

inline CCPi::cgls_base::cgls_base(const int niterations)
  : iterations(niterations)
{
}

inline CCPi::tv_regularization::tv_regularization(const real alph,
						  const real t, const real L,
						  const real mu, const int c)
  : alpha(alph), tau(t), init_L(L), init_mu(mu), constraint(c)
{
}

#endif // CCPI_RECON_ALGORITHMS
