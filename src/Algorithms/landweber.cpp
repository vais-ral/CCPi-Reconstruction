
#include <iostream>
#include "base_types.hpp"
#include "instruments.hpp"
#include "algorithms.hpp"
#include "timer.hpp"
#include "ui_calls.hpp"
#include "blas.hpp"
#include "landweber.hpp"

#ifndef USE_TIMER
#  define USE_TIMER false
#endif // USE_TIMER

bool CCPi::landweberLS::reconstruct(instrument *device, voxel_data &voxels,
				    const real origin[3],
				    const real voxel_size[3])
{
  const voxel_data::size_type *sz = voxels.shape();
  //sl_int n_vox = sl_int(sz[0]) * sl_int(sz[1]) * sl_int(sz[2]);
  //voxel_type *const x = voxels.data();
  pixel_data &b = device->get_pixel_data();
  int n_angles = device->get_num_angles();
  int n_h = device->get_num_h_pixels();
  int n_v = device->get_num_v_pixels();
  sl_int nx = sl_int(sz[0]);
  sl_int ny = sl_int(sz[1]);
  sl_int nz = sl_int(sz[2]);

  initialise_progress(2 * iterations, "LandweberLS iterating...");

  // Iterate.
  timer iter_time(USE_TIMER);
  for (int iter = 0; iter < iterations; iter++) {
    iter_time.reset();
    pixel_data r(boost::extents[n_angles][n_h][n_v]);
    // r = A * x;
    device->forward_project(r, voxels, origin, voxel_size,
			    (int)sz[0], (int)sz[1], (int)sz[2]);
    // r = b - r; -> r = b - A * x
    scal_xby(b, -1.0, r, n_angles, n_h, n_v);
    update_progress(2 * iter + 1);
    voxel_data s(boost::extents[sz[0]][sz[1]][sz[2]],
		 boost::c_storage_order());
    // s = At * r
    device->backward_project(r, s, origin, voxel_size,
			     (int)sz[0], (int)sz[1], (int)sz[2]);
    // x = x + lam * s; -> x = x + lam * (At * r)
    sum_axpy(lambda, s, voxels, nx, ny, nz);
    // project onto positive subspace
    clamp_min(voxels, 0.0, nx, ny, nz);
    update_progress(2 * iter + 2);
    iter_time.accumulate();
    iter_time.output("Iteration ");
  }
  return true;
}

bool CCPi::landweberLS::supports_blocks() const
{
  return false;
}
