
#ifndef CCPI_RECON_TVREG
#define CCPI_RECON_TVREG

namespace CCPi {

  class tv_regularization : public reconstruction_alg {
  public:
    tv_regularization(const real alph, const real t, const real L,
		      const real mu, const int c);

    bool reconstruct(instrument *device, voxel_data &voxels,
		     const real origin[3], const real voxel_size[3]);
    bool supports_blocks() const;

    static void tvreg_core(voxel_data &xkp1, real &fxkp1, real &hxkp1,
			   real &gxkp1, real_1dr &fxkp1l, int &kend,
			   const real voxel_size[], const pixel_data &b,
			   const real alpha, real tau, real bL, real bmu,
			   real epsb_rel,int k_max, const int Ddim,
			   const int Dm, const int Dn, const int Dl,
			   const sl_int prodDims, const int ctype,
			   real_1dr &d, real_1dr &c, const bool ghxl,
			   const bool xl, real_1dr &hxkp1l, real_1dr &gxkp1l,
			   real_1dr &xlist, const bool verbose,
			   int &numGrad, int &numBack, int &numFunc,
			   int &numRest, real_1dr &Lklist, real_1dr &muklist,
			   std::list<int> &rp, const real grid_offset[],
			   class instrument *device);

  private:
    real alpha;
    real tau;
    real init_L;
    real init_mu;
    int constraint;

    bool reconstruct(class instrument *device, pixel_data &b,
		     voxel_data &voxels, const real origin[3],
		     const real voxel_size[3]);
  };

}

inline CCPi::tv_regularization::tv_regularization(const real alph,
						  const real t, const real L,
						  const real mu, const int c)
  : alpha(alph), tau(t), init_L(L), init_mu(mu), constraint(c)
{
}

#endif // CCPI_RECON_TVREG
