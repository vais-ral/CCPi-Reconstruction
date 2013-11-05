
#include <iostream>

#include "base_types.hpp"
#include "mpi.hpp"
#include "utils.hpp"
#include "instruments.hpp"
#include "algorithms.hpp"
#include "results.hpp"

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
  bool fast_projection = false;
  int niterations = 5;
  CCPi::devices device = CCPi::dev_Nikon_XTek;
  //CCPi::devices device = CCPi::dev_Diamond_I12;
  CCPi::algorithms algorithm = CCPi::alg_CGLS;
  CCPi::instrument *instrument = 0;
  CCPi::output_format write_format = CCPi::bgs_float_dump;
  bool clamp_output = true;
  std::string output_name = "phantom";
  std::string data_file =
    "/home/bgs/scratch/ccpi/Bird_skull/Bird_skull_2001.xtekct";
  const int pixels_per_voxel = 1;
  // vertical size to break data up into for processing
  const int blocking_factor = 250;
  // number of GPUs etc if using accelerated code
  //const int num_devices = 1;
  // Todo - improve for TVReg
  const real alpha = 0.005;
  const real tau = 1e-4;
  const real l = 6930;
  const real mu = 0.5;
  const real tv_reg_constraint = 3;
  // distributed memory
  machine::initialise();
  int num_processors = machine::get_number_of_processors();
  // Todo - get stuff rather than the above test defaults here
  switch (device) {
  case CCPi::dev_Diamond_I12:
    instrument = new CCPi::Diamond;
    break;
  case CCPi::dev_Nikon_XTek:
    instrument = new CCPi::Nikon_XTek;
    break;
  default:
    std::cerr << "ERROR: Unknown device type\n";
    break;
  }
  if (num_processors == 1 or instrument->supports_distributed_memory()) {
    std::string path;
    std::string filename;
    CCPi::split_path_and_name(data_file, path, filename);
    if (instrument->setup_experimental_geometry(path, filename, phantom)) {
      if (instrument->read_data_size(path, phantom)) {
	// calculate blocks
	if (instrument->get_num_h_pixels() % pixels_per_voxel != 0)
	  std::cerr << "Number of horizontal pixels doesn't match voxels "
		    << instrument->get_num_h_pixels() << '\n';
	int nx_voxels = instrument->get_num_h_pixels() / pixels_per_voxel;
	int ny_voxels = nx_voxels;
	if (instrument->get_num_v_pixels() % pixels_per_voxel != 0)
	  std::cerr << "Number of vertical pixels doesn't match voxels "
		    << instrument->get_num_v_pixels() << '\n';
	int maxz_voxels = instrument->get_num_v_pixels() / pixels_per_voxel;
	int nz_voxels = 0;
	int block_size = 0;
	int block_step = 0;
	if (blocking_factor == 0 and num_processors == 1) {
	  nz_voxels = instrument->get_num_v_pixels() / pixels_per_voxel;
	  block_size = nz_voxels;
	  block_step = nz_voxels;
	} else if (instrument->supports_blocks()) {
	  int sz = 1;
	  if (blocking_factor > 0)
	    sz = blocking_factor;
	  int n_vox = instrument->get_num_v_pixels() / pixels_per_voxel;
	  if (n_vox / (sz * num_processors) < 1)
	    std::cerr << "Reduce blocking factor or number of processors\n";
	  block_size = sz;
	  block_step = block_size * num_processors;
	  nz_voxels = block_size;
	} else if (num_processors == 1) {
	  std::cerr << "Ignoring blocking factor - not supported by device\n";
	  nz_voxels = instrument->get_num_v_pixels() / pixels_per_voxel;
	  block_size = nz_voxels;
	  block_step = nz_voxels;
	}
	int z_data_size = block_size * pixels_per_voxel;
	int z_data_step = block_step * pixels_per_voxel;
	instrument->set_v_block(z_data_size);
	int block_offset = machine::get_processor_id() * block_size;
	int z_data_offset = block_offset * pixels_per_voxel;
	voxel_data voxels(boost::extents[nx_voxels][ny_voxels][nz_voxels],
			  boost::fortran_storage_order());
	real full_vox_origin[3];
	real voxel_size[3];
	if (instrument->finish_voxel_geometry(full_vox_origin, voxel_size,
					      nx_voxels, ny_voxels,
					      maxz_voxels)) {
	  // can modify offsets and end if parallel beam to solve subregion
	  // xTodo - if v_pixels % pix_pe_vox != 0 will this cause problems?
	  int end_value = instrument->total_num_v_pixels();
	  bool ok = false;
	  bool first = true;
	  do {
	    ok = false;
	    if (block_offset + block_size > nz_voxels)
	      block_size = nz_voxels - block_offset;
	    if (z_data_offset + z_data_size > end_value) {
	      z_data_size = end_value - z_data_offset;
	      instrument->set_v_block(z_data_size);
	    }
	    real voxel_origin[3];
	    voxel_origin[0] = full_vox_origin[0];
	    voxel_origin[1] = full_vox_origin[1];
	    voxel_origin[2] = full_vox_origin[2] + block_offset * voxel_size[2];
	    if (instrument->read_scans(path, z_data_offset,
				       z_data_size, first, phantom)) {
	      for (int i = 0; i < nz_voxels; i++) {
		for (int j = 0; j < ny_voxels; j++) {
		  for (int k = 0; k < nx_voxels; k++) {
		    voxels[k][j][i] = 0.0;
		  }
		}
	      }
	      if (beam_harden)
		instrument->apply_beam_hardening();
	      if (fast_projection and first)
		instrument->setup_projection_matrix(voxel_origin, voxel_size,
						    nx_voxels, ny_voxels,
						    nz_voxels);
	      switch (algorithm) {
	      case CCPi::alg_FDK:
		std::cerr << "ERROR: FDK not implmented - Todo\n";
		break;
	      case CCPi::alg_CGLS:
		ok = CCPi::cgls_reconstruction(instrument, voxels, voxel_origin,
					       voxel_size, niterations);
		break;
	      case CCPi::alg_TVreg:
		ok = CCPi::tv_regularization(instrument, voxels, voxel_origin,
					     voxel_size, alpha, tau, l, mu,
					     tv_reg_constraint);
		break;
	      default:
		std::cerr << "ERROR: Unknown algorithm\n";
		break;
	      }
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
    }
  } else
    std::cerr 
      << "Program does not support distributed memory for this instrument\n";
  machine::exit();
}
