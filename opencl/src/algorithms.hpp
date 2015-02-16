
#ifndef CCPI_RECON_ALGORITHMS
#define CCPI_RECON_ALGORITHMS

#include <list>

namespace CCPi {

  enum algorithms { alg_FDK, alg_CGLS, alg_TVreg, alg_BiCGLS, alg_BiCGSTABLS };

  class reconstruction_alg {
  public:
    virtual ~reconstruction_alg();
    virtual bool reconstruct(class instrument *device, voxel_data &voxels,
			     const real origin[3],
			     const real voxel_size[3]) = 0;
    virtual bool supports_blocks() const = 0;
  };

  class cgls_base : public reconstruction_alg {
  public:
    cgls_base(const int niterations);

    bool reconstruct(class instrument *device, voxel_data &voxels,
		     const real origin[3], const real voxel_size[3]);

  protected:
    virtual int size_of_voxel_norm(const int nz) const = 0;
    virtual void normalise_voxels(voxel_data &v, const sl_int nx,
				  const sl_int ny, const sl_int nz,
				  voxel_1d &norm) const = 0;
    virtual void pixel_update(const pixel_data &Ad, pixel_data &b,
			      const sl_int n_angles, const sl_int n_v,
			      const sl_int n_h, const voxel_data &d, 
			      voxel_data &v, const sl_int nx,
			      const sl_int ny, const sl_int nz,
			      const voxel_1d &norm) const = 0;
    virtual void voxel_update(const voxel_data &s, voxel_data &v,
			      const sl_int nx, const sl_int ny, const sl_int nz,
			      voxel_1d &norm) const = 0;
    int get_iterations() const;

  private:
    int iterations;
  };

  class cgls_3d : public cgls_base {
  public:
    cgls_3d(const int niterations);

    bool supports_blocks() const;

  protected:
    int size_of_voxel_norm(const int nz) const;
    void normalise_voxels(voxel_data &v, const sl_int nx, const sl_int ny,
			  const sl_int nz, voxel_1d &norm) const;
    void pixel_update(const pixel_data &Ad, pixel_data &b,
		      const sl_int n_angles, const sl_int n_v,
		      const sl_int n_h, const voxel_data &d, voxel_data &voxels,
		      const sl_int nx, const sl_int ny, const sl_int nz,
		      const voxel_1d &norm) const;
    void voxel_update(const voxel_data &s, voxel_data &v, const sl_int nx,
		      const sl_int ny, const sl_int nz, voxel_1d &norm) const;
  };

  class cgls_2d : public cgls_base {
  public:
    cgls_2d(const int niterations, const int ppv);

    bool supports_blocks() const;

  protected:
    int size_of_voxel_norm(const int nz) const;
    void normalise_voxels(voxel_data &v, const sl_int nx, const sl_int ny,
			  const sl_int nz, voxel_1d &norm) const;
    void pixel_update(const pixel_data &Ad, pixel_data &b,
		      const sl_int n_angles, const sl_int n_v,
		      const sl_int n_h, const voxel_data &d, voxel_data &voxels,
		      const sl_int nx, const sl_int ny, const sl_int nz,
		      const voxel_1d &norm) const;
    void voxel_update(const voxel_data &s, voxel_data &v, const sl_int nx,
		      const sl_int ny, const sl_int nz, voxel_1d &norm) const;

  private:
    int pixels_per_voxel;
  };

  class bi_cgls_3d : public cgls_3d {
  public:
    bi_cgls_3d(const int niterations);

    bool reconstruct(class instrument *device, voxel_data &voxels,
		     const real origin[3], const real voxel_size[3]);

  private:
    voxel_type voxel_update(voxel_data &v, voxel_data &p0,
			    voxel_data &r0, const voxel_data &vq,
			    const sl_int nx, const sl_int ny, const sl_int nz,
			    const pixel_data &q, const sl_int n_angles,
			    const sl_int n_h, const sl_int n_v,
			    const voxel_type gamma) const;
  };

  class bi_cgstabls_3d : public cgls_3d {
  public:
    bi_cgstabls_3d(const int niterations);

    bool reconstruct(class instrument *device, voxel_data &voxels,
		     const real origin[3], const real voxel_size[3]);
  };

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

inline CCPi::reconstruction_alg::~reconstruction_alg()
{
}

inline CCPi::cgls_base::cgls_base(const int niterations)
  : iterations(niterations)
{
}

inline int CCPi::cgls_base::get_iterations() const
{
  return iterations;
}

inline CCPi::cgls_3d::cgls_3d(const int niterations)
  : cgls_base(niterations)
{
}

inline CCPi::cgls_2d::cgls_2d(const int niterations, const int ppv)
  : cgls_base(niterations), pixels_per_voxel(ppv)
{
}

inline CCPi::tv_regularization::tv_regularization(const real alph,
						  const real t, const real L,
						  const real mu, const int c)
  : alpha(alph), tau(t), init_L(L), init_mu(mu), constraint(c)
{
}

inline CCPi::bi_cgls_3d::bi_cgls_3d(const int niterations)
  : cgls_3d(niterations)
{
}

inline CCPi::bi_cgstabls_3d::bi_cgstabls_3d(const int niterations)
  : cgls_3d(niterations)
{
}

#endif // CCPI_RECON_ALGORITHMS
