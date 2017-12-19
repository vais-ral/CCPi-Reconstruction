
#ifndef CCPI_RECON_INSTRUMENTS
#define CCPI_RECON_INSTRUMENTS

#include "CCPiDefines.h"
#include <map>

namespace CCPi {

  enum devices { dev_Diamond_I12, dev_Nikon_XTek };

  enum ring_artefact_alg { no_ring_artefacts, ring_artefacts_column,
			   ring_artefacts_aml };

  void CCPI_EXPORT remove_column_ring_artefacts(pixel_data &pixels, const sl_int nangles,
				    const sl_int nh, const sl_int nv);
  void CCPI_EXPORT remove_aml_ring_artefacts(pixel_data &pixels, const sl_int nangles,
				 const sl_int nh, const sl_int nv,
				 const real param_n, const real param_r,
				 const int num_series);

  class CCPI_EXPORT instrument {
  public:
    instrument();
    virtual ~instrument();

    virtual bool setup_experimental_geometry(const std::string path,
					     const std::string file,
					     const real rotation_centre,
					     const int pixels_per_voxel,
					     const bool phantom = false) = 0;
    virtual bool setup_experimental_geometry(const numpy_3d &pix_array,
					     const numpy_1d &angle_array,
					     const real rotation_centre,
					     const int pixels_per_voxel) = 0;
    virtual bool setup_experimental_geometry(const numpy_3d &pix_array,
					     const numpy_1d &angle_array,
					     const numpy_1d &h_offsets,
					     const numpy_1d &v_offsets,
					     const int pixels_per_voxel,
					     const real source_x,
					     const real detector_x,
					     const real pixel_h_size,
					     const real pixel_v_size,
					     const real m_radius,
					     const bool has_offsets) = 0;
    virtual bool read_scans(const std::string path, const int offset,
			    const int block_size, const bool first,
			    const bool phantom = false) = 0;
    virtual bool read_scans(const numpy_3d &pixel_array, const int offset,
			    const int block_size, bool is_pixel_in_log) = 0;
    int get_num_h_pixels() const;
    int get_num_v_pixels() const;
    int total_num_v_pixels() const;
    int get_num_angles() const;
    virtual bool finish_voxel_geometry(real voxel_origin[3], real voxel_size[3],
				       const int nx, const int ny,
				       const int nz) const = 0;
    virtual void get_xy_size(int &nx, int &ny, const int pixels_per_voxel) = 0;
    virtual void apply_beam_hardening() = 0;
    virtual void forward_project(pixel_data &pixels, voxel_data &voxels,
				 const real origin[3], const real width[3],
				 const int nx, const int ny, const int nz) = 0;
    virtual void backward_project(pixel_data &pixels, voxel_data &voxels,
				  const real origin[3], const real width[3],
				  const int nx, const int ny, const int nz) = 0;
    virtual void backward_project(voxel_data &voxels,
				  const real origin[3], const real width[3],
				  const int nx, const int ny, const int nz) = 0;

    virtual bool supports_distributed_memory() const = 0;
    virtual bool supports_blocks() const = 0;
    void set_v_block(const int size);

    // parallel beam
    static void forward_project(const real_1d &det_y, const real_1d &det_z,
				const real_1d &phi,
				pixel_data &ray_data, voxel_data &vol_data,
				const int n_angles, const int n_rays_y,
				const int n_rays_z, const real grid_offset[3],
				const real voxel_size[3], const int nx_voxels,
				const int ny_voxels, const int nz_voxels);
    static void backward_project(const real_1d &det_y, const real_1d &det_z,
				 const real_1d &phi,
				 pixel_data &ray_data, voxel_data &vol_data,
				 const int n_angles, const int n_rays_y,
				 const int n_rays_z, const real grid_offset[3],
				 const real voxel_size[3], const int nx_voxels,
				 const int ny_voxels, const int nz_voxels);
    // cone beam
    static void forward_project(const real source_x, const real source_y,
				const real source_z, const real det_x,
				const real_1d &det_y, const real_1d &det_z,
				const real_1d &phi,
				pixel_data &ray_data, voxel_data &vol_data,
				const int n_angles, const int n_rays_y,
				const int n_rays_z, const real grid_offset[3],
				const real voxel_size[3], const int nx_voxels,
				const int ny_voxels, const int nz_voxels);
    static void backward_project(const real source_x, const real source_y,
				 const real source_z, const real det_x,
				 const real_1d &det_y, const real_1d &det_z,
				 const real_1d &phi,
				 pixel_data &ray_data, voxel_data &vol_data,
				 const int n_angles, const int n_rays_y,
				 const int n_rays_z, const real grid_offset[3],
				 const real voxel_size[3], const int nx_voxels,
				 const int ny_voxels, const int nz_voxels);

    // Todo - protect these? CGLS uses them
    pixel_3d &get_pixel_data();
    const real_1d &get_phi() const; //expose this to read angles
	
  protected:
    pixel_3d &create_pixel_data();

    const real_1d &get_h_pixels() const;
    const real_1d &get_v_pixels() const;
    const real_1d &get_all_v_pixels() const;
    int get_data_v_size() const;
    int get_data_v_offset() const;

    real_1d &set_phi(const int n);
    real_1d &set_h_pixels(const int n);
    real_1d &set_v_pixels(const int n);
    void set_v_offset(const int offset);
    void adjust_h_pixels(const real centre);
    int calc_v_alignment(const int n, const int pix_per_vox, const bool cone);

  private:
    real_1d phi;
    real_1d horizontal_pixels;
    real_1d vertical_pixels;
    real_1d all_vertical_pixels;
    int n_angles;
    int n_horizontal_pixels;
    int n_vertical_pixels;
    int total_vertical_pixels;
    int v_offset;
    int data_v_size;
    int data_v_offset;
    pixel_3d *pixels;
  };

  class CCPI_EXPORT cone_beam : public instrument {
  public:
    void forward_project(pixel_data &pixels, voxel_data &voxels,
			 const real origin[3], const real width[3],
			 const int nx, const int ny, const int nz);
    void backward_project(pixel_data &pixels, voxel_data &voxels,
			  const real origin[3], const real width[3],
			  const int nx, const int ny, const int nz);
    void backward_project(voxel_data &voxels,
			  const real origin[3], const real width[3],
			  const int nx, const int ny, const int nz);

    bool supports_distributed_memory() const;
    bool supports_blocks() const;

    // Kludge for Matlab interface.
    void set_params(const real sx, const real sy, const real sz, const real dx,
		    real dy[], real dz[], real ang[], const int ny,
		    const int nz, const int nang);
    static void f2D(const real source_x, const real source_y,
		    const real source_z, const real detector_x,
		    const real_1d &h_pixels, const real_1d &v_pixels,
		    const real_1d &angles, pixel_data &pixels,
		    voxel_data &voxels, const int n_angles, const int n_h,
		    const int n_v, const real grid_offset[3],
		    const real voxel_size[3], const int nx_voxels,
		    const int ny_voxels, const int nz_voxels);
    static void b2D(const real source_x, const real source_y,
		    const real source_z, const real detector_x,
		    const real_1d &h_pixels, const real_1d &v_pixels,
		    const real_1d &angles, pixel_data &pixels,
		    voxel_data &voxels, const int n_angles, const int n_h,
		    const int n_v, const real vox_origin[3],
		    const real vox_size[3], const int nx, const int ny,
		    const int nz, const bool limited_memory = false);

    real get_source_x() const;
    real get_source_y() const;
    real get_source_z() const;
    real get_detector_x() const;

  protected:
    void set_source(const real x, const real y, const real z);
    void set_detector(const real x);

    void safe_forward_project(pixel_data &pixels, voxel_data &voxels,
			      const real origin[3], const real width[3],
			      const int nx, const int ny, const int nz);

  private:
    real source_x;
    real source_y;
    real source_z;
    // Todo - does this need further generalisation?
    real detector_x;

    static void calc_xy_z(pixel_type *const pixels, const voxel_ptr_1d &voxels,
			  const recon_1d &alpha_xy, const int n,
			  const recon_type pzbz, const recon_type inv_dz,
			  const int nv, const int nz, const int midp,
			  const recon_1d &delta_z, const recon_1d &inv_delz,
			  const recon_1d &vox_z);
    static void calc_ah_z(const pixel_ptr_1d &pixels, voxel_type *const voxels,
			  const recon_1d &alpha_xy_0,
			  const recon_1d &alpha_xy_1, const int n,
			  const recon_type pzbz, const recon_type inv_dz,
			  const int nv, const int nz, const int midp,
			  const recon_1d &delta_z, const recon_1d &inv_delz,
			  const recon_1d &vox_z, const recon_type pzdv,
			  const recon_type z_1, const recon_type z_nm);
    static void fproject_xy(const real p1_x, const real p1_y, const real p2_x,
			    const real p2_y, pixel_data &pixels,
			    voxel_data &voxels, const real b_x, const real b_y,
			    const real d_x, const real d_y, const int nx,
			    const int ny, const int nz, const int a,
			    const int h, const int nv, const int midp,
			    const recon_1d &delta_z, const recon_1d &inv_delz,
			    const recon_1d &vox_z, const recon_type pzbz,
			    const recon_type inv_dz, const sl_int ij_base,
			    const sl_int nyz);
    static void bproject_ah(const real source_x, const real source_y,
			    pixel_data &pixels, voxel_data &voxels,
			    const real x_0, const real y_0, const real x_n,
			    const real y_n, const int nz, const int i,
			    const int j, const int n_angles, const int n_h,
			    const int n_v, const real_1d &h_pixels,
			    const int midp, const real_1d &cangle,
			    const real_1d &sangle, const recon_1d &delta_z,
			    const recon_1d &inv_delz, const recon_1d &vox_z,
			    const recon_type pzbz, const recon_type inv_dz,
			    const recon_type pzdv, const recon_type z_1,
			    const recon_type z_nm, const real_1d &p1x,
			    const real_1d &p1y, const real_1d &cdetx,
			    const real_1d &sdetx, const real_1d &ilcphi,
			    const real_1d &ilsphi, const int a_off);
    static void b2D(const real source_x, const real source_y,
		    const real source_z, const real detector_x,
		    const real_1d &h_pixels, const real_1d &v_pixels,
		    const real_1d &angles, pixel_data &pixels,
		    voxel_data &voxels, const int n_angles, const int n_h,
		    const int n_v, const real vox_origin[3],
		    const real vox_size[3], const int nx, const int ny,
		    const int nz, const recon_2d &d_conv);
  };

  class CCPI_EXPORT parallel_beam : public instrument {
  public:
    void forward_project(pixel_data &pixels, voxel_data &voxels,
			 const real origin[3], const real width[3],
			 const int nx, const int ny, const int nz);
    void backward_project(pixel_data &pixels, voxel_data &voxels,
			  const real origin[3], const real width[3],
			  const int nx, const int ny, const int nz);
    void backward_project(voxel_data &voxels,
			  const real origin[3], const real width[3],
			  const int nx, const int ny, const int nz);

    bool supports_distributed_memory() const;
    bool supports_blocks() const;

    // for matlab interface
    static void f2D(const real_1d &h_pixels, const real_1d &v_pixels,
		    const real_1d &angles, const int n_angles,
		    const int nh_pixels, const int nv_pixels,
		    const real vox_origin[3], const real vox_size[3],
		    const int nx, const int ny, const int nz,
		    pixel_data &pixels, voxel_data &voxels);
    static void b2D(const real_1d &h_pixels, const real_1d &v_pixels,
		    const real_1d &angles, pixel_data &pixels,
		    voxel_data &voxels, const int n_angles, const int nh_pixels,
		    const int nv_pixels, const real grid_offset[3],
		    const real voxel_size[3], const int nx_voxels,
		    const int ny_voxels, const int nz_voxels);

  protected:
    void safe_forward_project(pixel_data &pixels, voxel_data &voxels,
			      const real origin[3], const real width[3],
			      const int nx, const int ny, const int nz);

  private:
    static void calc_xy_z(pixel_type *const pixels, const voxel_ptr_1d &voxels,
			  const recon_1d &l_xy, const int n, const int nv,
			  const int nz, const int_1d &mapping,
			  const int map_type);
    static void calc_ah_z(const pixel_ptr_1d &pixels, voxel_type *const voxels,
			  const recon_1d &l_xy, const int n, const int nv,
			  const int nz, const int_1d &mapping,
			  const int map_type);
    static void fproject_xy(const real p2_x, const real p2_y,
			    pixel_data &pixels, voxel_data &voxels,
			    const real b_x, const real b_y, const real d_x,
			    const real d_y, const int nx, const int ny,
			    const int nz, const int a, const int h,
			    const int nv, const recon_type d_conv,
			    const real cphi, const real sphi,
			    const sl_int ij_base, const sl_int nyz,
			    const int_1d &mapping, const int map_type);
    static void bproject_ah(pixel_data &pixels, voxel_data &voxels,
			    const real x_0, const real y_0, const real x_n,
			    const real y_n, const real d_x, const real d_y,
			    const int nz, const int i, const int j,
			    const int n_angles, const int n_h, const int n_v,
			    const real_1d &h_pixels, const real_1d &cangle,
			    const real_1d &sangle, const real_1d &y_offset,
			    const real_1d &i_offset, const real_1d &length,
			    const real h_pix0, const real ihp_step,
			    const int a_off, const int_1d &mapping,
			    const int map_type);
    static void gen_mapping(int_1d &mapping, int &map_type,
			    const real_1d &v_pixels, const real vox_z,
			    const real size_z, const int nv);
  };

  class CCPI_EXPORT Diamond : public parallel_beam {
  public:
    Diamond(const bool use_hp = false, const real hpj = 0.0,
	    const int hp_np = 1, const bool use_ring = false,
	    const real ra_n = 0.0, const real ra_r = 0.0, const int ra_ns = 0,
	    const ring_artefact_alg ra_alg = no_ring_artefacts);
    bool setup_experimental_geometry(const std::string path,
				     const std::string file,
				     const real rotation_centre,
				     const int pixels_per_voxel,
				     const bool phantom);
    bool setup_experimental_geometry(const numpy_3d &pix_array,
				     const numpy_1d &angle_array,
				     const real rotation_centre,
				     const int pixels_per_voxel);
    bool setup_experimental_geometry(const numpy_3d &pix_array,
				     const numpy_1d &angle_array,
				     const numpy_1d &h_offsets,
				     const numpy_1d &v_offsets,
				     const int pixels_per_voxel,
				     const real source_x, const real detector_x,
				     const real pixel_h_size,
				     const real pixel_v_size,
				     const real m_radius,
				     const bool has_offsets);
    bool read_scans(const std::string path, const int offset,
		    const int block_size, const bool first, const bool phantom);
    bool read_scans(const numpy_3d &pixel_array, const int offset,
		    const int block_size,bool is_pixel_in_log);
    bool finish_voxel_geometry(real voxel_origin[3], real voxel_size[3],
			       const int nx, const int ny, const int nz) const;
    void get_xy_size(int &nx, int &ny, const int pixels_per_voxel);
    void apply_beam_hardening();

  private:
    real hp_jump;
    int hp_num_pix;
    real aml_param_n;
    real aml_param_r;
    int ra_num_series;
    ring_artefact_alg ra_algorithm;
    bool use_high_peaks;
    bool use_ring_artefacts;
    std::string name;
    real h_vox_size;
    real v_vox_size;

    bool create_phantom();
    bool build_phantom(const int offset, const int block_size);
    bool read_data_size(const std::string path, const real rotation_centre,
			const int pixels_per_voxel);
    bool read_data(const std::string path, const int offset,
		   const int block_size, const bool first);
    void high_peaks_before(const real jump, const int num_pix);
    void ring_artefact_removal(const ring_artefact_alg alg, const real param_n,
			       const real param_r, const int num_series);
  };

  class CCPI_EXPORT Nikon_XTek : public cone_beam {
  public:
    bool setup_experimental_geometry(const std::string path,
				     const std::string file,
				     const real rotation_centre,
				     const int pixels_per_voxel,
				     const bool phantom);
    bool setup_experimental_geometry(const numpy_3d &pix_array,
				     const numpy_1d &angle_array,
				     const real rotation_centre,
				     const int pixels_per_voxel);
    bool setup_experimental_geometry(const numpy_3d &pix_array,
				     const numpy_1d &angle_array,
				     const numpy_1d &h_offsets,
				     const numpy_1d &v_offsets,
				     const int pixels_per_voxel,
				     const real source_x, const real detector_x,
				     const real pixel_h_size,
				     const real pixel_v_size,
				     const real m_radius,
				     const bool has_offsets);
    bool read_scans(const std::string path, const int offset,
		    const int block_size, const bool first, const bool phantom);
    bool read_scans(const numpy_3d &pixel_array, const int offset,
		    const int block_size,bool is_pixel_in_log);
    bool finish_voxel_geometry(real voxel_origin[3], real voxel_size[3],
			       const int nx, const int ny, const int nz) const;
    void get_xy_size(int &nx, int &ny, const int pixels_per_voxel);
	real get_mask_radius(){return mask_radius;}
    void apply_beam_hardening();
	void initialise_phantom();
    void initialise_withconfig(const std::string& config_path);
    voxel_data initialise_phantom_volume(real_1d &h_pixels, real_1d &v_pixels);
	
  private:
    real mask_radius;
    real white_level;
    std::string basename;
    real h_vox_size;
    real v_vox_size;
    real scattering;
    real coeff_x4;
    real coeff_x3;
    real coeff_x2;
    real coeff_x1;
    real coeff_x0;
    real scale;

    bool create_phantom();
    bool build_phantom();
    bool read_config_file(const std::string path, const std::string file,
			  const int pixels_per_voxel);
    bool read_angles(const std::string path, const real init_angle,
		     const int n);
    bool read_images(const std::string path);
    bool find_centre(const int v_slice);
  };

}

inline CCPi::instrument::instrument() : pixels(0)
{
}

inline const real_1d &CCPi::instrument::get_phi() const
{
  return phi;
}

inline const real_1d &CCPi::instrument::get_h_pixels() const
{
  return horizontal_pixels;
}

inline const real_1d &CCPi::instrument::get_v_pixels() const
{
  if (v_offset == 0)
    return all_vertical_pixels;
  else
    return vertical_pixels;
}

inline const real_1d &CCPi::instrument::get_all_v_pixels() const
{
  return all_vertical_pixels;
}

inline real_1d &CCPi::instrument::set_phi(const int n)
{
  n_angles = n;
  phi.resize(n);
  return phi;
}

inline real_1d &CCPi::instrument::set_h_pixels(const int n)
{
  n_horizontal_pixels = n;
  horizontal_pixels.resize(n);
  return horizontal_pixels;
}

inline real_1d &CCPi::instrument::set_v_pixels(const int n)
{
  n_vertical_pixels = n;
  total_vertical_pixels = n;
  v_offset = 0;
  vertical_pixels.resize(n);
  all_vertical_pixels.resize(n);
  return all_vertical_pixels;
}

inline int CCPi::instrument::get_num_angles() const
{
  return n_angles;
}

inline int CCPi::instrument::get_num_h_pixels() const
{
  return n_horizontal_pixels;
}

inline int CCPi::instrument::get_num_v_pixels() const
{
  return n_vertical_pixels;
}

inline int CCPi::instrument::total_num_v_pixels() const
{
  return total_vertical_pixels;
}

inline void CCPi::instrument::adjust_h_pixels(const real centre)
{
  for (int i = 0; i < n_horizontal_pixels; i++)
    horizontal_pixels[i] += centre;
}

inline void CCPi::instrument::set_v_block(const int size)
{
  n_vertical_pixels = size;
}

inline void CCPi::instrument::set_v_offset(const int offset)
{
  v_offset = offset;
  for (int i = 0; i < total_vertical_pixels - offset; i++)
    vertical_pixels[i] = all_vertical_pixels[i + offset];
}

inline int CCPi::instrument::get_data_v_size() const
{
  return data_v_size;
}

inline int CCPi::instrument::get_data_v_offset() const
{
  return data_v_offset;
}

inline real CCPi::cone_beam::get_source_x() const
{
  return source_x;
}

inline real CCPi::cone_beam::get_source_y() const
{
  return source_y;
}

inline real CCPi::cone_beam::get_source_z() const
{
  return source_z;
}

inline real CCPi::cone_beam::get_detector_x() const
{
  return detector_x;
}

inline void CCPi::cone_beam::set_source(const real x, const real y,
					const real z)
{
  source_x = x;
  source_y = y;
  source_z = z;
}

inline void CCPi::cone_beam::set_detector(const real x)
{
  detector_x = x;
}

inline CCPi::Diamond::Diamond(const bool use_hp, const real hpj,
			      const int hp_np, const bool use_ring,
			      const real ra_n, const real ra_r, const int ra_ns,
			      const ring_artefact_alg ra_alg)
  : hp_jump(hpj), hp_num_pix(hp_np), aml_param_n(ra_n), aml_param_r(ra_r),
    ra_num_series(ra_ns), ra_algorithm(ra_alg), use_high_peaks(use_hp),
    use_ring_artefacts(use_ring)
{
}

#endif // CCPI_RECON_INSTRUMENTS
