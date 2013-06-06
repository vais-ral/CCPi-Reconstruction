
#include <iostream>

#include "base_types.hpp"
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
  bool phantom = false;
  bool beam_harden = false;
  bool fast_projection = false;
  int niterations = 5;
  CCPi::devices device = CCPi::dev_Nikon_XTek;
  CCPi::algorithms algorithm = CCPi::alg_CGLS;
  CCPi::instrument *instrument = 0;
  CCPi::output_format write_format = CCPi::bgs_float_dump;
  bool clamp_output = true;
  std::string output_name = "phantom";
  std::string data_file =
    "/home/bgs/scratch/ccpi/Bird_skull/Bird_skull_2001.xtekct";
  const int pixels_per_voxel = 4;
  // Todo - improve for TVReg
  const real alpha = 0.005;
  const real tau = 1e-4;
  const real l = 6930;
  const real mu = 0.5;
  const real tv_reg_constraint = 3;
  // Todo - get stuff rather than the above test defaults here
  switch (device) {
  case CCPi::dev_Diamond_I13:
    instrument = new CCPi::Diamond;
    break;
  case CCPi::dev_Nikon_XTek:
    instrument = new CCPi::Nikon_XTek;
    break;
  default:
    std::cerr << "ERROR: Unknown device type\n";
    break;
  }
  std::string path;
  std::string filename;
  CCPi::split_path_and_name(data_file, path, filename);
  if (instrument->setup_experimental_geometry(path, filename, phantom)) {
    if (instrument->read_scans(path, phantom)) {
      if (instrument->get_num_h_pixels() % pixels_per_voxel != 0)
	std::cerr << "Number of horizontal pixels doesn't match voxels "
		  << instrument->get_num_h_pixels() << '\n';
      int nx_voxels = instrument->get_num_h_pixels() / pixels_per_voxel;
      int ny_voxels = nx_voxels;
      if (instrument->get_num_v_pixels() % pixels_per_voxel != 0)
	std::cerr << "Number of vertical pixels doesn't match voxels "
		  << instrument->get_num_v_pixels() << '\n';
      int nz_voxels = instrument->get_num_v_pixels() / pixels_per_voxel;
      voxel_data voxels(boost::extents[nx_voxels][ny_voxels][nz_voxels],
			boost::fortran_storage_order());
      for (int i = 0; i < nz_voxels; i++) {
	for (int j = 0; j < ny_voxels; j++) {
	  for (int k = 0; k < nx_voxels; k++) {
	    voxels[k][j][i] = 0.0;
	  }
	}
      }
      real voxel_origin[3];
      real voxel_size[3];
      if (instrument->finish_voxel_geometry(voxel_origin, voxel_size,
					    voxels)) {      
	if (beam_harden)
	  instrument->apply_beam_hardening();
	if (fast_projection)
	  instrument->setup_projection_matrix(voxel_origin, voxel_size,
					      nx_voxels, ny_voxels, nz_voxels);
	bool ok = false;
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
	  CCPi::write_results(output_name, voxels, voxel_origin,
			      voxel_size, write_format, clamp_output);
      }
    }
  }
}
