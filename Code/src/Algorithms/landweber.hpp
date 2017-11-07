
#ifndef CCPI_RECON_LANDWEBER
#define CCPI_RECON_LANDWEBER

namespace CCPi {

  class landweberLS : public reconstruction_alg {
  public:
    landweberLS(const int niterations, const real lam);

    bool reconstruct(class instrument *device, voxel_data &voxels,
		     const real origin[3], const real voxel_size[3]);
    bool supports_blocks() const;

  protected:
    int get_iterations() const;

  private:
    int iterations;
    real lambda;
  };

}

inline CCPi::landweberLS::landweberLS(const int niterations, const real lam)
  : iterations(niterations), lambda(lam)
{
}

inline int CCPi::landweberLS::get_iterations() const
{
  return iterations;
}

#endif // CCPI_RECON_LANDWEBER
