
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
	       boost::c_storage_order());
  init_data(d, nx, ny, nz);
  initialise_progress(2 * iterations + 1, "CGLS iterating...");
  device->backward_project(d, origin, voxel_size,
			   (int)sz[0], (int)sz[1], (int)sz[2]);

  voxel_1d normr2(size_of_voxel_norm(nz));
  normalise_voxels(d, nx, ny, nz, normr2);
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
      pixel_update(Ad, b, n_angles, n_v, n_h, d, voxels, nx, ny, nz, normr2);
    }
    update_progress(2 * iter + 2);
    {
      voxel_data s(boost::extents[sz[0]][sz[1]][sz[2]],
		   boost::c_storage_order());
      init_data(s, nx, ny, nz);
      device->backward_project(b, s, origin, voxel_size,
			       (int)sz[0], (int)sz[1], (int)sz[2]);

      // Update d vector.
      voxel_update(s, d, nx, ny, nz, normr2);
    }
    update_progress(2 * iter + 3);
    iter_time.accumulate();
    iter_time.output("Iteration ");
  }
  //delete [] d;
  return true;
}

int CCPi::cgls_3d::size_of_voxel_norm(const int nz) const
{
  return 1;
}

int CCPi::cgls_2d::size_of_voxel_norm(const int nz) const
{
  return nz;
}

void CCPi::cgls_3d::normalise_voxels(voxel_data &v, const sl_int nx,
				     const sl_int ny, const sl_int nz,
				     voxel_1d &norm) const
{
  norm[0] = norm_voxels(v, nx, ny, nz);
}

void CCPi::cgls_2d::normalise_voxels(voxel_data &v, const sl_int nx,
				     const sl_int ny, const sl_int nz,
				     voxel_1d &norm) const
{
  norm_voxels(v, nx, ny, nz, norm);
}

void CCPi::cgls_3d::pixel_update(const pixel_data &Ad, pixel_data &b,
				 const sl_int n_angles, const sl_int n_v,
				 const sl_int n_h, const voxel_data &d,
				 voxel_data &voxels, const sl_int nx,
				 const sl_int ny, const sl_int nz,
				 const voxel_1d &norm) const
{
  pixel_type alpha = norm_pixels(Ad, n_angles, n_v, n_h);
  alpha = norm[0] / alpha;
  sum_axpy(alpha, d, voxels, nx, ny, nz);
  sum_axpy(-alpha, Ad, b, n_angles, n_v, n_h);
}

void CCPi::cgls_2d::pixel_update(const pixel_data &Ad, pixel_data &b,
				 const sl_int n_angles, const sl_int n_v,
				 const sl_int n_h, const voxel_data &d,
				 voxel_data &voxels, const sl_int nx,
				 const sl_int ny, const sl_int nz,
				 const voxel_1d &norm) const
{
  pixel_1d alpha_v(n_v);
  norm_pixels(Ad, n_angles, n_v, n_h, alpha_v);
  pixel_1d alpha(nz);
  int count = 0;
  for (int i = 0; i < nz; i++) {
    int step = pixels_per_voxel;
    if (count + step > n_v)
      step =  n_v - count;
    alpha[i] = 0.0;
    for (int j = 0; j < step; j++)
      alpha[i] += alpha_v[count + j];
    count += step;
    alpha[i] = norm[i] / alpha[i];
  }
  sum_axpy(alpha, d, voxels, nx, ny, nz);
  sub_axpy(alpha, Ad, b, n_angles, n_v, n_h, pixels_per_voxel);
}

void CCPi::cgls_3d::voxel_update(const voxel_data &s, voxel_data &d,
				 const sl_int nx, const sl_int ny,
				 const sl_int nz, voxel_1d &norm) const
{
  real normr2_new = norm_voxels(s, nx, ny, nz);
  real beta = normr2_new / norm[0];
  norm[0] = normr2_new;
  scal_xby(s, beta, d, nx, ny, nz);
}

void CCPi::cgls_2d::voxel_update(const voxel_data &s, voxel_data &d,
				 const sl_int nx, const sl_int ny,
				 const sl_int nz, voxel_1d &norm) const
{
  voxel_1d normr2_new(nz);
  norm_voxels(s, nx, ny, nz, normr2_new);
  for (int i = 0; i < nz; i++) {
    voxel_type n = normr2_new[i];
    normr2_new[i] /= norm[i];
    norm[i] = n;
  }
  scal_xby(s, normr2_new, d, nx, ny, nz);
}

bool CCPi::cgls_3d::supports_blocks() const
{
  return false;
}

bool CCPi::cgls_2d::supports_blocks() const
{
  return true;
}
