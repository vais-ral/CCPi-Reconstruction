
#ifndef CCPI_RECON_INSTRUMENTS
#define CCPI_RECON_INSTRUMENTS

namespace CCPi {

  enum devices { dev_Diamond_I13, dev_Nikon_XTek };

  class instrument {
  public:
    virtual bool setup_experimental_geometry(const std::string path,
					     const std::string file,
					     const bool phantom = false) = 0;
    virtual bool read_scans(const std::string path,
			    const bool phantom = false) = 0;
    virtual bool finish_voxel_geometry(real voxel_origin[3], real voxel_size[3],
				       const voxel_data &voxels) const = 0;
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

    template <class pixel_t, class voxel_t>
    static void forward_project(const real source_x, const real source_y,
				const real source_z, const real det_x,
				const real det_y[], const real det_z[],
				const real phi[], const real theta[],
				pixel_t ray_data[], voxel_t *const vol_data,
				const int n_angles, const int n_rays_y,
				const int n_rays_z, const real grid_offset[3],
				const real voxel_size[3], const int nx_voxels,
				const int ny_voxels, const int nz_voxels);
    template <class pixel_t, class voxel_t>
    static void backward_project(const real source_x, const real source_y,
				 const real source_z, const real det_x,
				 const real det_y[], const real det_z[],
				 const real phi[], const real theta[],
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

    real *get_phi() const;
    real *get_theta() const;
    real *get_h_pixels() const;
    real *get_v_pixels() const;
    int get_num_angles() const;
    int get_num_h_pixels() const;
    int get_num_v_pixels() const;

    void set_h_pixels(real *h_pixels, const int n);
    void set_v_pixels(real *v_pixels, const int n);
    void set_angles(real *p, real *t, const int n);
    void set_phi(real *p, const int n);

  private:
    real *phi;
    real *theta;
    real *horizontal_pixels;
    real *vertical_pixels;
    int n_angles;
    int n_horizontal_pixels;
    int n_vertical_pixels;
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
    // Todo - does this further generalisation?
    real detector_x;
  };

  class parallel_beam : public instrument {
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
  };

  class Diamond : public parallel_beam {
  public:
    bool setup_experimental_geometry(const std::string path,
				     const std::string file,
				     const bool phantom);
    bool read_scans(const std::string path, const bool phantom);
    bool finish_voxel_geometry(real voxel_origin[3], real voxel_size[3],
			       const voxel_data &voxels) const;
    void apply_beam_hardening();
  };

  class Nikon_XTek : public cone_beam {
  public:
    bool setup_experimental_geometry(const std::string path,
				     const std::string file,
				     const bool phantom);
    bool read_scans(const std::string path, const bool phantom);
    bool finish_voxel_geometry(real voxel_origin[3], real voxel_size[3],
			       const voxel_data &voxels) const;
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
  };

}

inline real *CCPi::instrument::get_phi() const
{
  return phi;
}

inline real *CCPi::instrument::get_theta() const
{
  return theta;
}

inline real *CCPi::instrument::get_h_pixels() const
{
  return horizontal_pixels;
}

inline real *CCPi::instrument::get_v_pixels() const
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

inline void CCPi::instrument::set_h_pixels(real *h_pixels, const int n)
{
  horizontal_pixels = h_pixels;
  n_horizontal_pixels = n;
}

inline void CCPi::instrument::set_v_pixels(real *v_pixels, const int n)
{
  vertical_pixels = v_pixels;
  n_vertical_pixels = n;
}

inline void CCPi::instrument::set_angles(real *p, real *t, const int n)
{
  phi = p;
  theta = t;
  n_angles = n;
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

#endif // CCPI_RECON_INSTRUMENTS
