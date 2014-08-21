
#include "base_types.hpp"
#include "utils.hpp"
#include "fbp.hpp"
#include "instruments.hpp"
#include "voxels.hpp"
#include "ui_calls.hpp"

void calculate_block_sizes(int &nx_voxels, int &ny_voxels, int &nz_voxels,
			   int &maxz_voxels, int &block_size, int &block_step,
			   const int num_processors, const int blocking_factor,
			   const int pixels_per_voxel,
			   CCPi::instrument *instrument,
			   const bool recon_blocks)
{
  // calculate blocks
  instrument->get_xy_size(nx_voxels, ny_voxels, pixels_per_voxel);
  maxz_voxels = instrument->get_num_v_pixels() / pixels_per_voxel;
  if (instrument->get_num_v_pixels() % pixels_per_voxel != 0)
    maxz_voxels++;
  nz_voxels = 0;
  block_size = 0;
  block_step = 0;
  if (blocking_factor == 0 and num_processors == 1) {
    nz_voxels = maxz_voxels;
    block_size = nz_voxels;
    block_step = nz_voxels;
  } else if (instrument->supports_blocks() and recon_blocks) {
    int sz = 1;
    if (blocking_factor > 0)
      sz = blocking_factor;
    if (maxz_voxels / (sz * num_processors) < 1)
      report_error("Reduce blocking factor or number of processors");
    block_size = sz;
    block_step = block_size * num_processors;
    nz_voxels = block_size;
  } else if (num_processors == 1) {
    report_error("Ignoring blocking factor - not supported by device");
    nz_voxels = maxz_voxels;
    block_size = nz_voxels;
    block_step = nz_voxels;
  }
}
