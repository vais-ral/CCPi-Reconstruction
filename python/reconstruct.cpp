
#include "base_types.hpp"
#include "mpi.hpp"
#include "utils.hpp"
#include "instruments.hpp"
#include "algorithms.hpp"
#include "results.hpp"
#include "voxels.hpp"

void reconstruct_cgls(std::string data_file, std::string output_base,
		      double rotation_centre, int resolution,
		      int niterations)
{
  // Todo ring artefacts choice etc.

  bool phantom = false;
  bool beam_harden = false;
  CCPi::output_format write_format = CCPi::native_dump;
  bool clamp_output = true;
  // vertical size to break data up into for processing
  const int blocking_factor = 0;
  // number of GPUs etc if using accelerated code
  //const int num_devices = 1;
  CCPi::instrument *instrument = new CCPi::Diamond(true, 0.05, 2, true, 0.000,
						   0.0050, 1,
						   CCPi::ring_artefacts_aml);
  CCPi::reconstruction_alg *algorithm = new CCPi::cgls_3d(niterations);
  //if (blocking_factor > 0 and instrument->supports_blocks())
  //  recon_algorithm = new CCPi::cgls_2d(niterations, pixels_per_voxel);
  machine::initialise();
  reconstruct(instrument, algorithm, data_file, output_base, rotation_centre,
	      resolution, blocking_factor, beam_harden, write_format,
	      clamp_output, phantom);
  machine::exit();
}

void reconstruct_tvreg()
{
}
