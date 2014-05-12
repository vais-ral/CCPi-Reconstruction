
#include <iostream>
#include "base_types.hpp"
#include "instruments.hpp"
#include "algorithms.hpp"
#include "timer.hpp"
#include "ui_calls.hpp"
#include "blas.hpp"

#ifndef USE_TIMER
#  define USE_TIMER false
#endif // USE_TIMER

bool CCPi::cgls_base::reconstruct(instrument *device, voxel_data &voxels,
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

  // Prepare for CG iteration.
  voxel_data d(boost::extents[sz[0]][sz[1]][sz[2]],
	       boost::fortran_storage_order());
  init_data(d, nx, ny, nz);
  initialise_progress(2 * iterations + 1, "CGLS iterating...");
  device->backward_project(d, origin, voxel_size,
			   (int)sz[0], (int)sz[1], (int)sz[2]);

  real normr2 = norm_voxels(d, nx, ny, nz);
  update_progress(1);

  // Iterate.
  timer iter_time(USE_TIMER);
  for (int iter = 0; iter < iterations; iter++) {
    //add_output("iter ");
    //add_output(j + 1);
    //send_output();
    iter_time.reset();
    // Update x and r vectors.
    {
      pixel_data Ad(boost::extents[n_angles][n_v][n_h]);
      init_data(Ad, n_angles, n_v, n_h);
      device->forward_project(Ad, d, origin, voxel_size,
			      (int)sz[0], (int)sz[1], (int)sz[2]);
      real alpha = norm_pixels(Ad, n_angles, n_v, n_h);
      alpha = normr2 / alpha;
      sum_axpy(alpha, d, voxels, nx, ny, nz);
      sum_axpy(-alpha, Ad, b, n_angles, n_v, n_h);
    }
    update_progress(2 * iter + 2);
    {
      voxel_data s(boost::extents[sz[0]][sz[1]][sz[2]],
		   boost::fortran_storage_order());
      init_data(s, nx, ny, nz);
      device->backward_project(b, s, origin, voxel_size,
			       (int)sz[0], (int)sz[1], (int)sz[2]);

      // Update d vector.
      real normr2_new = norm_voxels(s, nx, ny, nz);
      real beta = normr2_new / normr2;
      normr2 = normr2_new;
      scal_xby(s, beta, d, nx, ny, nz);
    }
    update_progress(2 * iter + 3);
    iter_time.accumulate();
    iter_time.output("Iteration ");
  }
  //delete [] d;
  return true;
}
