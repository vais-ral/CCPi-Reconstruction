
#ifndef CCPI_RECON_SIRT
#define CCPI_RECON_SIRT

namespace CCPi {

  class sirt : public reconstruction_alg {
  public:
    sirt(const int niterations);

    bool reconstruct(class instrument *device, voxel_data &voxels,
		     const real origin[3], const real voxel_size[3]);
    bool supports_blocks() const;

  protected:
    int get_iterations() const;

  private:
    int iterations;
  };

}

inline CCPi::sirt::sirt(const int niterations)
  : iterations(niterations)
{
}

inline int CCPi::sirt::get_iterations() const
{
  return iterations;
}

#endif // CCPI_RECON_SIRT
