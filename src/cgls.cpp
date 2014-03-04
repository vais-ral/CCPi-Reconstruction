
#include <iostream>
#include "base_types.hpp"
#include "fbp.hpp"
#include "instruments.hpp"
#include "algorithms.hpp"
#include "timer.hpp"

#ifndef USE_TIMER
#  define USE_TIMER false
#endif // USE_TIMER

bool CCPi::cgls_reconstruction(const instrument *device, voxel_data &voxels,
			       const real origin[3], const real voxel_size[3],
			       const int iterations)
{
  const voxel_data::size_type *sz = voxels.shape();
  sl_int n_vox = sl_int(sz[0]) * sl_int(sz[1]) * sl_int(sz[2]);
  voxel_type *const x = voxels.data();
  pixel_type *const b = device->get_pixel_data();

  // Prepare for CG iteration.
  std::cout << "Preparing for CG iteration...\n";
  voxel_type *d = new voxel_type[n_vox];
  for (sl_int i = 0; i < n_vox; i++)
    d[i] = 0.0;
  device->backward_project(d, origin, voxel_size,
			   (int)sz[0], (int)sz[1], (int)sz[2]);
  sl_int n_rays = device->get_data_size();

  real normr2 = 0.0;
  for (sl_int i = 0; i < n_vox; i++)
    normr2 += d[i] * d[i];

  // work space
  void *temp = 0;
  std::size_t psize = sizeof(pixel_type);
  std::size_t vsize = sizeof(voxel_type);  
  if (n_rays * psize > n_vox * vsize)
    temp = new pixel_type[n_rays];
  else
    temp = new voxel_type[n_vox];
  pixel_type *Ad = (pixel_type *)temp;
  voxel_type *s = (voxel_type *)temp;

  // Iterate.
  timer iter_time(USE_TIMER);
  for (int j = 0; j < iterations; j++) {
    std::cout << "iter " << j + 1 << '\n';
    iter_time.reset();
    // Update x and r vectors.
    for (sl_int i = 0; i < n_rays; i++)
      Ad[i] = 0.0;
    device->forward_project(Ad, d, origin, voxel_size,
			    (int)sz[0], (int)sz[1], (int)sz[2]);
    real alpha = 0.0;
    for (sl_int i = 0; i < n_rays; i++)
      alpha += Ad[i] * Ad[i];
    alpha = normr2 / alpha;
    for (sl_int i = 0; i < n_vox; i++)
      x[i] += alpha * d[i];
    for (sl_int i = 0; i < n_rays; i++)
      b[i] -= alpha * Ad[i];
    for (sl_int i = 0; i < n_vox; i++)
      s[i] = 0.0;
    device->backward_project(b, s, origin, voxel_size,
			     (int)sz[0], (int)sz[1], (int)sz[2]);

    // Update d vector.
    real normr2_new = 0.0;
    for (sl_int i = 0; i < n_vox; i++)
      normr2_new += s[i] * s[i];
    real beta = normr2_new / normr2;
    normr2 = normr2_new;
    for (sl_int i = 0; i < n_vox; i++)
      d[i] = s[i] + beta * d[i];
    iter_time.accumulate();
    iter_time.output("Iteration ");
  }
  delete [] s;
  delete [] d;
  return true;
}
