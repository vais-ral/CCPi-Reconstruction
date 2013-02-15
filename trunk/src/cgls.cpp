
#include <iostream>
#include "base_types.hpp"
#include "instruments.hpp"
#include "algorithms.hpp"

bool CCPi::cgls_reconstruction(const instrument *device, voxel_data &voxels,
			       const real origin[3], const real voxel_size[3],
			       const int iterations)
{
  const voxel_data::size_type *sz = voxels.shape();
  long n_vox = long(sz[0]) * long(sz[1]) * long(sz[2]);
  voxel_type *const x = voxels.data();
  pixel_type *const b = device->get_pixel_data();

  // Prepare for CG iteration.
  std::cout << "Preparing for CG iteration...\n";
  voxel_type *d = new voxel_type[n_vox];
  for (long i = 0; i < n_vox; i++)
    d[i] = 0.0;
  device->backward_project(d, origin, voxel_size,
			   (int)sz[0], (int)sz[1], (int)sz[2]);
  long n_rays = device->get_data_size();

  real normr2 = 0.0;
  for (long i = 0; i < n_vox; i++)
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
  for (int j = 0; j < iterations; j++) {
    std::cout << "iter " << j + 1 << '\n';
    // Update x and r vectors.
    for (long i = 0; i < n_rays; i++)
      Ad[i] = 0.0;
    device->forward_project(Ad, d, origin, voxel_size,
			    (int)sz[0], (int)sz[1], (int)sz[2]);
    real alpha = 0.0;
    for (long i = 0; i < n_rays; i++)
      alpha += Ad[i] * Ad[i];
    alpha = normr2 / alpha;
    for (long i = 0; i < n_vox; i++)
      x[i] += alpha * d[i];
    for (long i = 0; i < n_rays; i++)
      b[i] -= alpha * Ad[i];
    for (long i = 0; i < n_vox; i++)
      s[i] = 0.0;
    device->backward_project(b, s, origin, voxel_size,
			     (int)sz[0], (int)sz[1], (int)sz[2]);

    // Update d vector.
    real normr2_new = 0.0;
    for (long i = 0; i < n_vox; i++)
      normr2_new += s[i] * s[i];
    real beta = normr2_new / normr2;
    normr2 = normr2_new;
    for (long i = 0; i < n_vox; i++)
      d[i] = s[i] + beta * d[i];
  }
  delete [] s;
  delete [] d;
  return true;
}
