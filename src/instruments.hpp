
#ifndef CCPI_RECON_INSTRUMENTS
#define CCPI_RECON_INSTRUMENTS

namespace CCPi {

  enum devices { dev_Diamond_I13, dev_Nikon_XTek };

  class instrument {
  public:
    virtual bool setup_experimental_geometry(const std::string file,
					     const bool phantom = false) = 0;
    virtual bool read_scans(const bool phantom = false) = 0;
    virtual void forward_project(voxel_data &voxels, const real origin[3],
				 const real width[3]) = 0;
    virtual void backward_project(voxel_data &voxels, const real origin[3],
				  const real width[3]) = 0;
  };

  // Todo - do these serve any purpose?
  class cone_beam : public instrument{
  };

  class parallel_beam : public instrument {
  };

  class Diamond : public parallel_beam {
  public:
    bool setup_experimental_geometry(const std::string file,
				     const bool phantom);
    bool read_scans(const bool phantom);
    void forward_project(voxel_data &voxels, const real origin[3],
			 const real width[3]);
    void backward_project(voxel_data &voxels, const real origin[3],
			  const real width[3]);
  };

  class Nikon_XTek : public cone_beam {
  public:
    bool setup_experimental_geometry(const std::string file,
				     const bool phantom);
    bool read_scans(const bool phantom);
    void forward_project(voxel_data &voxels, const real origin[3],
			 const real width[3]);
    void backward_project(voxel_data &voxels, const real origin[3],
			  const real width[3]);

    template <class pixel_t, class voxel_t>
    static void forward_project(const real source_x, const real source_y,
				const real source_z, const real det_x,
				const real det_y[], const real det_z[],
				const real angles[], pixel_t ray_data[],
				voxel_t *const vol_data, const int n_angles,
				const int n_rays_y, const int n_rays_z,
				const real grid_offset[3],
				const real voxel_size[3], const int nx_voxels,
				const int ny_voxels, const int nz_voxels);
    template <class pixel_t, class voxel_t>
    static void backward_project(const real source_x, const real source_y,
				 const real source_z, const real det_x,
				 const real det_y[], const real det_z[],
				 const real angles[], pixel_t ray_data[],
				 voxel_t *const vol_data, const int n_angles,
				 const int n_rays_y, const int n_rays_z,
				 const real grid_offset[3],
				 const real voxel_size[3], const int nx_voxels,
				 const int ny_voxels, const int nz_voxels);

  private:
    real source_x;
    real source_y;
    real source_z;
    real detector_x;
    real *angles;
    real *horizontal_pixels;
    real *vertical_pixels;
    int n_angles;
    int n_horizontal_pixels;
    int n_vertical_pixels;
    pixel_type *pixels;

    bool create_phantom();
    bool build_phantom();
    bool read_config_file(const std::string file);
    bool read_images();
  };

}

#endif // CCPI_RECON_INSTRUMENTS
