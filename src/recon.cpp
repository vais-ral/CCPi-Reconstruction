
#include <iostream>

#include "base_types.hpp"
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
    - output format and location
    - ?
  */
  // Todo - usage messages if started up wrong?
  bool phantom = true;
  CCPi::devices device = CCPi::dev_Nikon_XTek;
  CCPi::algorithms algorithm = CCPi::alg_CGLS;
  CCPi::instrument *instrument = 0;
  std::string data_file;
  int nx_voxels = 500;
  int ny_voxels = 500;
  int nz_voxels = 500;
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
  voxel_data voxels(boost::extents[nx_voxels][ny_voxels][nz_voxels],
		    boost::fortran_storage_order());
  real voxel_origin[3];
  real voxel_size[3];
  std::cerr << "Todo - origin/voxel size\n";
  if (instrument->setup_experimental_geometry(data_file, phantom)) {
    if (instrument->read_scans(phantom)) {
      bool ok = false;
      switch (algorithm) {
      case CCPi::alg_FDK:
	std::cerr << "ERROR: FDK not implmented - Todo\n";
	break;
      case CCPi::alg_CGLS:
	ok = CCPi::cgls_reconstruction(instrument, voxels, voxel_origin,
				       voxel_size);
	break;
      case CCPi::alg_TVreg:
	std::cerr << "ERROR: TVreg not implmented - Todo\n";
	break;
      default:
	std::cerr << "ERROR: Unknown algorithm\n";
	break;
      }
      if (ok)
	CCPi::write_results(voxels);
    }
  }
}
