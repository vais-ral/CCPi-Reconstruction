
#ifndef CCPI_RECON_CGLS
#define CCPI_RECON_CGLS
#include "CCPiDefines.h"

namespace CCPi {

  class CCPI_EXPORT cgls_base : public reconstruction_alg {
  public:
    cgls_base(const int niterations);

    bool reconstruct(class instrument *device, voxel_data &voxels,
		     const real origin[3], const real voxel_size[3]);
    void convergence_data(real_1d &data) const;

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
    void set_norm(const real norm, const int iter);

  private:
    int iterations;
    real_1d norm_r;
  };

  class CCPI_EXPORT cgls_3d : public cgls_base {
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

  class CCPI_EXPORT cgls_2d : public cgls_base {
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

  class CCPI_EXPORT bi_cgls_3d : public cgls_3d {
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

  class CCPI_EXPORT bi_cgstabls_3d : public cgls_3d {
  public:
    bi_cgstabls_3d(const int niterations);

    bool reconstruct(class instrument *device, voxel_data &voxels,
		     const real origin[3], const real voxel_size[3]);
  };

  class CCPI_EXPORT cgls_regularize : public cgls_3d {
  public:
    cgls_regularize(const int niterations, const real param);
    bool reconstruct(class instrument *device, voxel_data &voxels,
		     const real origin[3], const real voxel_size[3]);

  private:
    real regularisation_param;

    virtual void regularize(voxel_data &b, const voxel_data &a,
			    const int nx, const int ny, const int nz) = 0;
  };

  class CCPI_EXPORT cgls_tikhonov : public cgls_regularize {
  public:
    cgls_tikhonov(const int niterations, const real param);

  private:
    void regularize(voxel_data &b, const voxel_data &a,
		    const int nx, const int ny, const int nz);
  };

  class CCPI_EXPORT cgls_tv_reg : public cgls_regularize {
  public:
    cgls_tv_reg(const int niterations, const real param);

  private:
    void regularize(voxel_data &b, const voxel_data &a,
		    const int nx, const int ny, const int nz);
  };

}

inline CCPi::cgls_base::cgls_base(const int niterations)
  : iterations(niterations), norm_r(iterations)
{
}

inline int CCPi::cgls_base::get_iterations() const
{
  return iterations;
}

inline void CCPi::cgls_base::set_norm(const real norm, const int iter)
{
  norm_r[iter] = norm;
}

inline CCPi::cgls_3d::cgls_3d(const int niterations)
  : cgls_base(niterations)
{
}

inline CCPi::cgls_2d::cgls_2d(const int niterations, const int ppv)
  : cgls_base(niterations), pixels_per_voxel(ppv)
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

inline CCPi::cgls_regularize::cgls_regularize(const int niterations,
					      const real param)
  : cgls_3d(niterations), regularisation_param(param)
{
}

inline CCPi::cgls_tikhonov::cgls_tikhonov(const int niterations,
					  const real param)
  : cgls_regularize(niterations, param)
{
}

inline CCPi::cgls_tv_reg::cgls_tv_reg(const int niterations,
				      const real param)
  : cgls_regularize(niterations, param)
{
}

#endif // CCPI_RECON_CGLS
