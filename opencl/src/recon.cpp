
#include <iostream>

#include "base_types.hpp"
#include "mpi.hpp"
#include "utils.hpp"
#include "instruments.hpp"
#include "algorithms.hpp"
#include "results.hpp"
#include "voxels.hpp"
#include "cgls.hpp"
#include "tv_reg.hpp"
#include "landweber.hpp"
#include "mlem.hpp"
#include "sirt.hpp"

int main()
{
  /*
    need input on
    - instrument type
    - algorithm to use
    - file/directory for data
    - size of reconstruction (number of voxels)
    - beam hardening
    - output format and location
    - ?
  */
  // Todo - usage messages if started up wrong?
  bool phantom = true;
  bool beam_harden = false;
  int niterations = 1;
  //CCPi::devices device = CCPi::dev_Nikon_XTek;
  CCPi::devices device = CCPi::dev_Diamond_I12;
  CCPi::algorithms algorithm = CCPi::alg_CGLS;
  CCPi::instrument *instrument = 0;
  CCPi::reconstruction_alg *recon_algorithm = 0;
  CCPi::output_format write_format = CCPi::bgs_float_dump;
  bool clamp_output = true;
  std::string output_name = "phantom";
  std::string data_file =
    "/home/bgs/scratch/ccpi/Bird_skull/Bird_skull_2001.xtekct";
  real rotation_centre = -1.0;
  const int pixels_per_voxel = 1;
  // vertical size to break data up into for processing
  const int blocking_factor = 0;
  // number of GPUs etc if using accelerated code
  //const int num_devices = 1;
  // Todo - improve for TVReg
  const real alpha = 0.005;
  const real tau = 1e-4;
  const real l = 6930;
  const real mu = 0.5;
  const real tv_reg_constraint = 3;
  // Todo - what should this be?
  const real lambda = 0.1;
  bool setup_ok = true;
  // Todo - get stuff rather than the above test defaults here
  switch (device) {
  case CCPi::dev_Diamond_I12:
    instrument = new CCPi::Diamond(true, 0.05, 2, true, 0.000, 0.0050, 1,
				   CCPi::ring_artefacts_aml);
    break;
  case CCPi::dev_Nikon_XTek:
    instrument = new CCPi::Nikon_XTek;
    break;
  default:
    std::cerr << "ERROR: Unknown device type\n";
    setup_ok = false;
    break;
  }
  switch (algorithm) {
  case CCPi::alg_FDK:
    std::cerr << "ERROR: FDK not implmented - Todo\n";
    setup_ok = false;
    break;
  case CCPi::alg_CGLS:
    if (blocking_factor > 0 and instrument->supports_blocks())
      recon_algorithm = new CCPi::cgls_2d(niterations, pixels_per_voxel);
    else
      recon_algorithm = new CCPi::cgls_3d(niterations);
    break;
  case CCPi::alg_landweber:
    recon_algorithm = new CCPi::landweberLS(niterations, lambda);
    break;
  case CCPi::alg_MLEM:
    recon_algorithm = new CCPi::mlem(niterations);
    break;
  case CCPi::alg_SIRT:
    recon_algorithm = new CCPi::sirt(niterations);
    break;
  case CCPi::alg_TVreg:
    recon_algorithm = new CCPi::tv_regularization(alpha, tau, l, mu,
						  tv_reg_constraint);
    break;
  case CCPi::alg_BiCGLS:
    recon_algorithm = new CCPi::bi_cgls_3d(niterations);
    break;
  case CCPi::alg_BiCGSTABLS:
    recon_algorithm = new CCPi::bi_cgstabls_3d(niterations);
    break;
  case CCPi::alg_CGLS_Tikhonov:
    recon_algorithm = new CCPi::cgls_tikhonov(niterations, 0.01);
    break;
  default:
    std::cerr << "ERROR: Unknown algorithm\n";
    setup_ok = false;
  }
  if (setup_ok) {
    // distributed memory
    machine::initialise();
    int num_processors = machine::get_number_of_processors();
    if (num_processors == 1 or instrument->supports_distributed_memory()) {
      real vox_origin[3];
      real vox_size[3];
      voxel_data *voxels = reconstruct(instrument, recon_algorithm, data_file,
				       output_name, vox_origin, vox_size,
				       rotation_centre, pixels_per_voxel,
				       blocking_factor, beam_harden,
				       write_format, clamp_output, phantom);
      delete voxels;
    } else
      std::cerr 
	<< "Program does not support distributed memory for this instrument\n";
    machine::exit();
  }
}
