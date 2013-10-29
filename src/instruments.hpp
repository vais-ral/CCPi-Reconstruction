
#ifndef CCPI_RECON_INSTRUMENTS
#define CCPI_RECON_INSTRUMENTS

#include <map>

namespace CCPi {

  struct map_index {
    int x;
    int y;

    bool operator < (const map_index &b) const;
  };

  typedef std::map<map_index, real> projection_map;

  enum devices { dev_Diamond_I12, dev_Nikon_XTek };

  class instrument {
  public:
    virtual bool setup_experimental_geometry(const std::string path,
					     const std::string file,
					     const bool phantom = false) = 0;
    virtual bool read_scans(const std::string path, const int offset,
			    const int block_size, const bool first,
			    const bool phantom = false) = 0;
    virtual bool read_data_size(const std::string path,
				const bool phantom = false) = 0;
    int get_num_h_pixels() const;
    int get_num_v_pixels() const;
    int total_num_v_pixels() const;
    virtual bool finish_voxel_geometry(real voxel_origin[3], real voxel_size[3],
				       const int nx, const int ny,
				       const int nz) const = 0;
    virtual void apply_beam_hardening() = 0;
    virtual void forward_project(pixel_type *pixels, voxel_type *const voxels,
				 const real origin[3], const real width[3],
				 const int nx, const int ny,
				 const int nz) const = 0;
    virtual void backward_project(pixel_type *pixels, voxel_type *const voxels,
				  const real origin[3], const real width[3],
				  const int nx, const int ny,
				  const int nz) const = 0;
    virtual void backward_project(voxel_type *const voxels,
				  const real origin[3], const real width[3],
				  const int nx, const int ny,
				  const int nz) const = 0;

    virtual void setup_projection_matrix(const real origin[3],
					 const real width[3],
					 const int nx, const int ny,
					 const int nz) = 0;

    virtual bool supports_distributed_memory() const = 0;
    virtual bool supports_blocks() const = 0;
    void set_v_block(const int size);

    // parallel beam
    template <class pixel_t, class voxel_t>
    static void forward_project(const real det_y[], const real det_z[],
				const real phi[],
				pixel_t ray_data[], voxel_t *const vol_data,
				const int n_angles, const int n_rays_y,
				const int n_rays_z, const real grid_offset[3],
				const real voxel_size[3], const int nx_voxels,
				const int ny_voxels, const int nz_voxels);
    template <class pixel_t, class voxel_t>
    static void backward_project(const real det_y[], const real det_z[],
				 const real phi[],
				 pixel_t ray_data[], voxel_t *const vol_data,
				 const int n_angles, const int n_rays_y,
				 const int n_rays_z, const real grid_offset[3],
				 const real voxel_size[3], const int nx_voxels,
				 const int ny_voxels, const int nz_voxels);
    // cone beam
    template <class pixel_t, class voxel_t>
    static void forward_project(const real source_x, const real source_y,
				const real source_z, const real det_x,
				const real det_y[], const real det_z[],
				const real phi[],
				pixel_t ray_data[], voxel_t *const vol_data,
				const int n_angles, const int n_rays_y,
				const int n_rays_z, const real grid_offset[3],
				const real voxel_size[3], const int nx_voxels,
				const int ny_voxels, const int nz_voxels);
    template <class pixel_t, class voxel_t>
    static void backward_project(const real source_x, const real source_y,
				 const real source_z, const real det_x,
				 const real det_y[], const real det_z[],
				 const real phi[],
				 pixel_t ray_data[], voxel_t *const vol_data,
				 const int n_angles, const int n_rays_y,
				 const int n_rays_z, const real grid_offset[3],
				 const real voxel_size[3], const int nx_voxels,
				 const int ny_voxels, const int nz_voxels);

    // Todo - protect these? CGLS uses them
    long get_data_size() const;
    pixel_type *const get_pixel_data() const;

  protected:
    pixel_type *create_pixel_data();
    void set_pixel_data(pixel_type *p, const long n);

    real *get_phi() const;
    real *get_h_pixels() const;
    real *get_v_pixels() const;
    real *get_all_v_pixels() const;
    int get_num_angles() const;

    void set_h_pixels(real *h_pixels, const int n);
    void set_v_pixels(real *v_pixels, const int n);
    void set_phi(real *p, const int n);
    void set_v_offset(const int offset);

  private:
    real *phi;
    real *horizontal_pixels;
    real *vertical_pixels;
    int n_angles;
    int n_horizontal_pixels;
    int n_vertical_pixels;
    int total_vertical_pixels;
    int v_offset;
    pixel_type *pixel_data;
  };

  class cone_beam : public instrument {
  public:
    void forward_project(pixel_type *pixels, voxel_type *const voxels,
			 const real origin[3], const real width[3],
			 const int nx, const int ny, const int nz) const;
    void backward_project(pixel_type *pixels, voxel_type *const voxels,
			  const real origin[3], const real width[3],
			  const int nx, const int ny, const int nz) const;
    void backward_project(voxel_type *const voxels,
			  const real origin[3], const real width[3],
			  const int nx, const int ny, const int nz) const;

    void setup_projection_matrix(const real origin[3], const real width[3],
				 const int nx, const int ny,
				 const int nz);

    bool supports_distributed_memory() const;
    bool supports_blocks() const;

    // Kludge for Matlab interface.
    void set_params(const real sx, const real sy, const real sz, const real dx,
		    real dy[], real dz[], real ang[], const int ny,
		    const int nz, const int nang);

  protected:
    real get_source_x() const;
    real get_source_y() const;
    real get_source_z() const;
    real get_detector_x() const;

    void set_source(const real x, const real y, const real z);
    void set_detector(const real x);

  private:
    real source_x;
    real source_y;
    real source_z;
    // Todo - does this need further generalisation?
    real detector_x;
  };

  class parallel_beam : public instrument {
  public:
    parallel_beam();

    void forward_project(pixel_type *pixels, voxel_type *const voxels,
			 const real origin[3], const real width[3],
			 const int nx, const int ny, const int nz) const;
    void backward_project(pixel_type *pixels, voxel_type *const voxels,
			  const real origin[3], const real width[3],
			  const int nx, const int ny, const int nz) const;
    void backward_project(voxel_type *const voxels,
			  const real origin[3], const real width[3],
			  const int nx, const int ny, const int nz) const;

    void setup_projection_matrix(const real origin[3], const real width[3],
				 const int nx, const int ny,
				 const int nz);

    bool supports_distributed_memory() const;
    bool supports_blocks() const;

  private:
    bool has_projection_matrix;
    long matrix_size;
    real *forward_matrix;
    long *forward_cols;
    long *forward_rows;
    real *backward_matrix;
    long *backward_cols;
    long *backward_rowb;
    long *backward_rowe;

    static void my_back_project(const real h_pixels[], const real v_pixels[],
				 const real angles[], pixel_type pixels[],
				 voxel_type *const voxels,
				 const int n_angles, const int nh_pixels,
				 const int nv_pixels, const real grid_offset[3],
				 const real voxel_size[3], const int nx_voxels,
				 const int ny_voxels, const int nz_voxels);
    static void map_2Dprojection(const real start[], const real end[],
				 const real b_x, const real b_y,
				 const real b_z, const real d_x,
				 const real d_y, const real d_z,
				 const int im_size_x, const int im_size_y,
				 const int im_size_z, const long z_offset,
				 projection_map &map);
    void setup_2D_matrix(const real det_y[], const real phi[],
			 const int n_angles, const int n_rays_z,
			 const int n_rays_y, const real grid_offset[3],
			 const real voxel_size[3], const int nx_voxels,
			 const int ny_voxels, const int nz_voxels);
    void forward_project_matrix(const real det_z[], pixel_type ray_data[],
				voxel_type *const vol_data, const int n_angles,
				const int n_rays_y, const int n_rays_z,
				const real grid_offset[3],
				const real voxel_size[3], const int nx_voxels,
				const int ny_voxels, const int nz_voxels) const;
    void backward_project_matrix(const real det_z[], pixel_type ray_data[],
				 voxel_type *const vol_data, const int n_angles,
				 const int n_rays_y, const int n_rays_z,
				 const real grid_offset[3],
				 const real voxel_size[3], const int nx_voxels,
				 const int ny_voxels,
				 const int nz_voxels) const;
  };

  class Diamond : public parallel_beam {
  public:
    bool setup_experimental_geometry(const std::string path,
				     const std::string file,
				     const bool phantom);
    bool read_scans(const std::string path, const int offset,
		    const int block_size, const bool first, const bool phantom);
    bool read_data_size(const std::string path, const bool phantom);
    bool finish_voxel_geometry(real voxel_origin[3], real voxel_size[3],
			       const int nx, const int ny, const int nz) const;
    void apply_beam_hardening();

  private:
    std::string name;

    bool create_phantom();
    bool build_phantom(const int offset, const int block_size);
    bool read_data(const std::string path, const int offset,
		   const int block_size, const bool first);
  };

  class Nikon_XTek : public cone_beam {
  public:
    bool setup_experimental_geometry(const std::string path,
				     const std::string file,
				     const bool phantom);
    bool read_scans(const std::string path, const int offset,
		    const int block_size, const bool first, const bool phantom);
    bool read_data_size(const std::string path, const bool phantom);
    bool finish_voxel_geometry(real voxel_origin[3], real voxel_size[3],
			       const int nx, const int ny, const int nz) const;
    void apply_beam_hardening();

  private:
    real offset[3];
    real mask_radius;
    real white_level;
    std::string basename;

    bool create_phantom();
    bool build_phantom();
    bool read_config_file(const std::string path, const std::string file);
    bool read_angles(const std::string datafile, const real init_angle,
		     const int n);
    bool read_images(const std::string path);
    void find_centre(const int v_slice);
  };

}

inline bool CCPi::map_index::operator < (const map_index &b) const
{
  if (x < b.x)
    return true;
  else if (x == b.x)
    return (y < b.y);
  else
    return false;
}

inline real *CCPi::instrument::get_phi() const
{
  return phi;
}

inline real *CCPi::instrument::get_h_pixels() const
{
  return horizontal_pixels;
}

inline real *CCPi::instrument::get_v_pixels() const
{
  return &vertical_pixels[v_offset];
}

inline real *CCPi::instrument::get_all_v_pixels() const
{
  return vertical_pixels;
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

inline void CCPi::instrument::set_h_pixels(real *h_pixels, const int n)
{
  horizontal_pixels = h_pixels;
  n_horizontal_pixels = n;
}

inline void CCPi::instrument::set_v_pixels(real *v_pixels, const int n)
{
  vertical_pixels = v_pixels;
  n_vertical_pixels = n;
  total_vertical_pixels = n;
  v_offset = 0;
}

inline void CCPi::instrument::set_v_block(const int size)
{
  n_vertical_pixels = size;
}

inline void CCPi::instrument::set_v_offset(const int offset)
{
  v_offset = offset;
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

inline CCPi::parallel_beam::parallel_beam()
  : has_projection_matrix(false), matrix_size(0), forward_matrix(0),
    forward_cols(0), forward_rows(0), backward_matrix(0), backward_cols(0),
    backward_rowb(0), backward_rowe(0)
{
}

#endif // CCPI_RECON_INSTRUMENTS
