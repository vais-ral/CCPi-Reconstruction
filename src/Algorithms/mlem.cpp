
#include <iostream>
#include "base_types.hpp"
#include "instruments.hpp"
#include "algorithms.hpp"
#include "timer.hpp"
#include "ui_calls.hpp"
#include "blas.hpp"
#include "mlem.hpp"

#ifndef USE_TIMER
#  define USE_TIMER false
#endif // USE_TIMER

bool CCPi::mlem::reconstruct(instrument *device, voxel_data &voxels,
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

  init_data(voxels, nx, ny, nz, voxel_type(1.0));

  // Prepare for iteration.
  initialise_progress(2 * iterations + 1, "MLEM iterating...");
  voxel_data normalisation(boost::extents[nx][ny][nz]);
  {
    pixel_data x_norm(boost::extents[n_angles][n_h][n_v]);
    init_data(x_norm, n_angles, n_h, n_v, pixel_type(1.0));
    device->backward_project(x_norm, normalisation, origin, voxel_size,
			     (int)sz[0], (int)sz[1], (int)sz[2]);
  }
  invert_x(normalisation, nx, ny, nz);
  pixel_data ratio(boost::extents[n_angles][n_h][n_v]);
  update_progress(1);

  // Iterate.
  timer iter_time(USE_TIMER);
  for (int iter = 0; iter < iterations; iter++) {
    iter_time.reset();
    {
      pixel_data sino_est(boost::extents[n_angles][n_h][n_v]);
      device->forward_project(sino_est, voxels, origin, voxel_size,
			      (int)sz[0], (int)sz[1], (int)sz[2]);
      // get ratio for values != 0 (> 0)?
      div_xyz(ratio, b, sino_est, n_angles, n_h, n_v);
    }
    /* Todo - the upper bound caused problems with the phantom. val <= 0.0?
    {
      pixel_1d val(n_v);
      for (int v = 0; v < n_v; v++)
	val[v] = 10.0 * ratio[n_angles / 2][n_h / 2][v]; // n/2+1?
      clampv_min_max(ratio, pixel_type(0.0), val, n_angles, n_h, n_v);
    }
    */
    // use this version of clamp instead of the commented out code above.
    clamp_min(ratio, 0.0, n_angles, n_h, n_v);
    update_progress(2 * iter + 2);
    {
      voxel_data bp_ratio(boost::extents[sz[0]][sz[1]][sz[2]],
			  boost::c_storage_order());
      device->backward_project(ratio, bp_ratio, origin, voxel_size,
			       (int)sz[0], (int)sz[1], (int)sz[2]);
      multsum_xyz(voxels, normalisation, bp_ratio, nx, ny, nz);
    }
    clamp_min(voxels, 0.0, nx, ny, nz);
    update_progress(2 * iter + 3);
    iter_time.accumulate();
    iter_time.output("Iteration ");
  }
  return true;
}

bool CCPi::mlem::supports_blocks() const
{
  return false;
}
