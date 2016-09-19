
#ifndef CCPI_RECON_MLEM
#define CCPI_RECON_MLEM

namespace CCPi {

  class mlem : public reconstruction_alg {
  public:
    mlem(const int niterations);

    bool reconstruct(class instrument *device, voxel_data &voxels,
		     const real origin[3], const real voxel_size[3]);
    bool supports_blocks() const;

  protected:
    int get_iterations() const;

  private:
    int iterations;
  };

}

inline CCPi::mlem::mlem(const int niterations)
  : iterations(niterations)
{
}

inline int CCPi::mlem::get_iterations() const
{
  return iterations;
}

#endif // CCPI_RECON_MLEM
