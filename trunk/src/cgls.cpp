
#include <iostream>
#include "base_types.hpp"
#include "instruments.hpp"
#include "algorithms.hpp"
#include "timer.hpp"
#include "ui_calls.hpp"

#ifndef USE_TIMER
#  define USE_TIMER false
#endif // USE_TIMER

bool CCPi::cgls_base::reconstruct(const instrument *device, voxel_data &voxels,
				  const real origin[3],
				  const real voxel_size[3])
{
  const voxel_data::size_type *sz = voxels.shape();
  //sl_int n_vox = sl_int(sz[0]) * sl_int(sz[1]) * sl_int(sz[2]);
  //voxel_type *const x = voxels.data();
  pixel_type *const b = device->get_pixel_data();

  // Prepare for CG iteration.
  voxel_data d(boost::extents[sz[0]][sz[1]][sz[2]],
	       boost::fortran_storage_order());
  //voxel_type *d = new voxel_type[n_vox];
  for (sl_int i = 0; i < sl_int(sz[2]); i++)
    for (sl_int j = 0; j < sl_int(sz[1]); j++)
      for (sl_int k = 0; k < sl_int(sz[0]); k++)
	d[i][j][k] = 0.0;
  initialise_progress(2 * iterations + 1, "CGLS iterating...");
  device->backward_project(d, origin, voxel_size,
			   (int)sz[0], (int)sz[1], (int)sz[2]);
  sl_int n_rays = device->get_data_size();

  real normr2 = 0.0;
  for (sl_int i = 0; i < sl_int(sz[2]); i++)
    for (sl_int j = 0; j < sl_int(sz[1]); j++)
      for (sl_int k = 0; k < sl_int(sz[0]); k++)
	normr2 += d[i][j][k] * d[i][j][k];
  update_progress(1);

  // Iterate.
  timer iter_time(USE_TIMER);
  for (int iter = 0; iter < iterations; iter++) {
    //add_output("iter ");
    //add_output(j + 1);
    //send_output();
    iter_time.reset();
    // Update x and r vectors.
    pixel_type *Ad = new pixel_type[n_rays];
    for (sl_int i = 0; i < n_rays; i++)
      Ad[i] = 0.0;
    device->forward_project(Ad, d, origin, voxel_size,
			    (int)sz[0], (int)sz[1], (int)sz[2]);
    real alpha = 0.0;
    for (sl_int i = 0; i < n_rays; i++)
      alpha += Ad[i] * Ad[i];
    alpha = normr2 / alpha;
    for (sl_int i = 0; i < sl_int(sz[2]); i++)
      for (sl_int j = 0; j < sl_int(sz[1]); j++)
	for (sl_int k = 0; k < sl_int(sz[0]); k++)
	  voxels[i][j][k] += alpha * d[i][j][k];
    for (sl_int i = 0; i < n_rays; i++)
      b[i] -= alpha * Ad[i];
    delete [] Ad;
    update_progress(2 * iter + 2);
    {
      voxel_data s(boost::extents[sz[0]][sz[1]][sz[2]],
		   boost::fortran_storage_order());
      for (sl_int i = 0; i < sl_int(sz[2]); i++)
	for (sl_int j = 0; j < sl_int(sz[1]); j++)
	  for (sl_int k = 0; k < sl_int(sz[0]); k++)
	    s[i][j][k] = 0.0;
      device->backward_project(b, s, origin, voxel_size,
			       (int)sz[0], (int)sz[1], (int)sz[2]);

      // Update d vector.
      real normr2_new = 0.0;
      for (sl_int i = 0; i < sl_int(sz[2]); i++)
	for (sl_int j = 0; j < sl_int(sz[1]); j++)
	  for (sl_int k = 0; k < sl_int(sz[0]); k++)
	    normr2_new += s[i][j][k] * s[i][j][k];
      real beta = normr2_new / normr2;
      normr2 = normr2_new;
      for (sl_int i = 0; i < sl_int(sz[2]); i++)
	for (sl_int j = 0; j < sl_int(sz[1]); j++)
	  for (sl_int k = 0; k < sl_int(sz[0]); k++)
	    d[i][j][k] = s[i][j][k] + beta * d[i][j][k];
    }
    update_progress(2 * iter + 3);
    iter_time.accumulate();
    iter_time.output("Iteration ");
  }
  //delete [] d;
  return true;
}
