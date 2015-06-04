
#ifndef CCPI_FBP
#define CCPI_FBP

namespace CCPi {

  class fbp_alg : public reconstruction_alg {
  public:
    fbp_alg(const filter_name_t label, const filter_window_t w,
	    const filter_norm_t n, const real b);

    bool reconstruct(class instrument *device, voxel_data &voxels,
		     const real origin[3], const real voxel_size[3]);
    bool supports_blocks() const;

  private:
    filter_name_t filter_name;
    filter_window_t window;
    filter_norm_t norm;
    real bandwidth;
  };

}

inline CCPi::fbp_alg::fbp_alg(const filter_name_t label,
			      const filter_window_t w,
			      const filter_norm_t n, const real b)
  : filter_name(label), window(w), norm(n), bandwidth(b)
{
}

#endif // CCPI_FBP
