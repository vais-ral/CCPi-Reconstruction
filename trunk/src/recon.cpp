
#include <iostream>

#include "base_types.hpp"
#include "mpi.hpp"
#include "utils.hpp"
#include "instruments.hpp"
#include "algorithms.hpp"
#include "results.hpp"
#include "voxels.hpp"
#include "blas.hpp"

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
  bool phantom = false;
  bool beam_harden = false;
  int niterations = 5;
  CCPi::devices device = CCPi::dev_Nikon_XTek;
  //CCPi::devices device = CCPi::dev_Diamond_I12;
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
  default:
    std::cerr << "ERROR: Unknown algorithm\n";
    setup_ok = false;
  }
  if (setup_ok) {
    // distributed memory
    machine::initialise();
    int num_processors = machine::get_number_of_processors();
    if (num_processors == 1 or instrument->supports_distributed_memory()) {
      std::string path;
      std::string filename;
      CCPi::split_path_and_name(data_file, path, filename);
      if (instrument->setup_experimental_geometry(path, filename,
						  rotation_centre,
						  pixels_per_voxel, phantom)) {
	int nx_voxels = 0;
	int ny_voxels = 0;
	int maxz_voxels = 0;
	int nz_voxels = 0;
	int block_size = 0;
	int block_step = 0;
	calculate_block_sizes(nx_voxels, ny_voxels, nz_voxels, maxz_voxels,
			      block_size, block_step, num_processors,
			      blocking_factor, pixels_per_voxel,
			      instrument, recon_algorithm->supports_blocks());
	int z_data_size = block_size * pixels_per_voxel;
	int z_data_step = block_step * pixels_per_voxel;
	instrument->set_v_block(z_data_size);
	int block_offset = machine::get_processor_id() * block_size;
	int z_data_offset = block_offset * pixels_per_voxel;
	real full_vox_origin[3];
	real voxel_size[3];
	if (instrument->finish_voxel_geometry(full_vox_origin, voxel_size,
					      nx_voxels, ny_voxels,
					      maxz_voxels)) {
	  // can modify offsets and end if parallel beam to solve subregion
	  int end_value = instrument->total_num_v_pixels();
	  bool ok = false;
	  bool first = true;
	  do {
	    ok = false;
	    if (block_offset + block_size > maxz_voxels)
	      block_size = maxz_voxels - block_offset;
	    if (z_data_offset + z_data_size > end_value) {
	      z_data_size = end_value - z_data_offset;
	      instrument->set_v_block(z_data_size);
	    }
	    nz_voxels = block_size;
	    real voxel_origin[3];
	    voxel_origin[0] = full_vox_origin[0];
	    voxel_origin[1] = full_vox_origin[1];
	    voxel_origin[2] = full_vox_origin[2]
	      + block_offset * voxel_size[2];
	    if (instrument->read_scans(path, z_data_offset,
				       z_data_size, first, phantom)) {
	      voxel_data voxels(boost::extents[nx_voxels][ny_voxels][nz_voxels]);
	      init_data(voxels, nx_voxels, ny_voxels, nz_voxels);
	      if (beam_harden)
		instrument->apply_beam_hardening();
	      ok = recon_algorithm->reconstruct(instrument, voxels,
						voxel_origin, voxel_size);
	      if (ok)
		CCPi::write_results(output_name, voxels, full_vox_origin,
				    voxel_size, block_offset, maxz_voxels,
				    write_format, clamp_output);
	    } else
	      ok = false;
	    first = false;
	    block_offset += block_step;
	    z_data_offset += z_data_step;
	  } while (ok and z_data_offset < end_value);
	}
      }
    } else
      std::cerr 
	<< "Program does not support distributed memory for this instrument\n";
    machine::exit();
  }
}
