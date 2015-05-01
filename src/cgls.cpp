
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
      pixel_data Ad(boost::extents[n_angles][n_h][n_v]);
      device->forward_project(Ad, d, origin, voxel_size,
			      (int)sz[0], (int)sz[1], (int)sz[2]);
      pixel_update(Ad, b, n_angles, n_v, n_h, d, voxels, nx, ny, nz, normr2);
    }
    update_progress(2 * iter + 2);
    {
      voxel_data s(boost::extents[sz[0]][sz[1]][sz[2]],
		   boost::c_storage_order());
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
  sum_axpy(-alpha, Ad, b, n_angles, n_h, n_v);
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

bool CCPi::bi_cgls_3d::reconstruct(instrument *device, voxel_data &voxels,
				   const real origin[3],
				   const real voxel_size[3])
{
  const voxel_data::size_type *sz = voxels.shape();
  //pixel_data &b = device->get_pixel_data();
  int n_angles = device->get_num_angles();
  int n_h = device->get_num_h_pixels();
  int n_v = device->get_num_v_pixels();
  sl_int nx = sl_int(sz[0]);
  sl_int ny = sl_int(sz[1]);
  sl_int nz = sl_int(sz[2]);

  // Prepare for CG iteration.
  voxel_data r0(boost::extents[sz[0]][sz[1]][sz[2]],
		boost::c_storage_order());
  initialise_progress(2 * get_iterations() + 1, "BiCGLS iterating...");
  device->backward_project(r0, origin, voxel_size,
			   (int)sz[0], (int)sz[1], (int)sz[2]);

  voxel_type gamma0 = 0.0;
  gamma0 = norm_voxels(r0, nx, ny, nz);
  //voxel_data rt0(boost::extents[sz[0]][sz[1]][sz[2]],
  //	 boost::c_storage_order());
  //copy(r0, rt0, nz, ny, nz);
  voxel_data p0(boost::extents[sz[0]][sz[1]][sz[2]],
		boost::c_storage_order());
  copy(r0, p0, nz, ny, nz);
  //voxel_data pt0(boost::extents[sz[0]][sz[1]][sz[2]],
  //	 boost::c_storage_order());  
  //copy(r0, pt0, nz, ny, nz);
  update_progress(1);

  // Issue - I don't see how pt0/rt0/qt can ever differ from p0/r0/q
  // Iterate.
  timer iter_time(USE_TIMER);
  for (int iter = 0; iter < get_iterations(); iter++) {
    iter_time.reset();
    pixel_data q(boost::extents[n_angles][n_h][n_v]);
    //pixel_data qt(boost::extents[n_angles][n_h][n_v]);
    //init_data(qt, n_angles, n_h, n_v);
    device->forward_project(q, p0, origin, voxel_size,
			    (int)sz[0], (int)sz[1], (int)sz[2]);
    //device->forward_project(qt, pt0, origin, voxel_size,
    //		      (int)sz[0], (int)sz[1], (int)sz[2]);
    update_progress(2 * iter + 2);
    voxel_data vq(boost::extents[sz[0]][sz[1]][sz[2]],
		  boost::c_storage_order());
    //voxel_data vqt(boost::extents[sz[0]][sz[1]][sz[2]],
    //	  boost::c_storage_order());
    //init_data(vqt, nx, ny, nz);
    device->backward_project(q, vq, origin, voxel_size,
			     (int)sz[0], (int)sz[1], (int)sz[2]);
    //device->backward_project(qt, vqt, origin, voxel_size,
    //		     (int)sz[0], (int)sz[1], (int)sz[2]);
    gamma0 = voxel_update(voxels, p0, r0, vq, nx, ny, nz,
			  q, n_angles, n_h, n_v, gamma0);
    update_progress(2 * iter + 3);
    iter_time.accumulate();
    iter_time.output("Iteration ");
  }
  return true;
}

voxel_type CCPi::bi_cgls_3d::voxel_update(voxel_data &v, voxel_data &p0,
					  voxel_data &r0, const voxel_data &vq,
					  const sl_int nx, const sl_int ny,
					  const sl_int nz, const pixel_data &q,
					  const sl_int n_angles,
					  const sl_int n_h, const sl_int n_v,
					  const voxel_type gamma) const
{
  real alpha = norm_pixels(q, n_angles, n_v, n_h);
  alpha = real(gamma) / alpha;
  // x = x + alpha * p0
  sum_axpy(alpha, p0, v, nx, ny, nz);
  // r0 = r0 - alpha * vq
  sum_axpy(- alpha, vq, r0, nx, ny, nz);
  real g = norm_voxels(r0, nx, ny, nz);
  real beta = g / real(gamma);
  scal_xby(p0, beta, r0, nx, ny, nz);
  return g;
}

bool CCPi::bi_cgstabls_3d::reconstruct(instrument *device, voxel_data &voxels,
				       const real origin[3],
				       const real voxel_size[3])
{
  const voxel_data::size_type *sz = voxels.shape();
  //pixel_data &b = device->get_pixel_data();
  int n_angles = device->get_num_angles();
  int n_h = device->get_num_h_pixels();
  int n_v = device->get_num_v_pixels();
  sl_int nx = sl_int(sz[0]);
  sl_int ny = sl_int(sz[1]);
  sl_int nz = sl_int(sz[2]);

  // Prepare for CG iteration.
  voxel_data r0(boost::extents[sz[0]][sz[1]][sz[2]],
		boost::c_storage_order());
  initialise_progress(2 * get_iterations() + 1, "BiCGSTABLS iterating...");
  device->backward_project(r0, origin, voxel_size,
			   (int)sz[0], (int)sz[1], (int)sz[2]);

  voxel_data r(boost::extents[sz[0]][sz[1]][sz[2]],
	       boost::c_storage_order());
  copy(r0, r, nx, ny, nz);
  voxel_data p0(boost::extents[sz[0]][sz[1]][sz[2]],
		boost::c_storage_order());
  copy(p0, r0, nx, ny, nz);
  voxel_data v0(boost::extents[sz[0]][sz[1]][sz[2]],
		boost::c_storage_order());
  copy(v0, r0, nx, ny, nz);
  voxel_type gamma0 = 1.0;
  voxel_type alpha = 1.0;
  voxel_type omega = 1.0;
  update_progress(1);

  // Iterate.
  timer iter_time(USE_TIMER);
  for (int iter = 0; iter < get_iterations(); iter++) {
    iter_time.reset();
    voxel_type gamma = norm_voxels(r0, r, nx, ny, nz);
    voxel_type beta = (gamma * alpha) / (gamma0 * omega);
    scal_xby(r, beta, p0, nx, ny, nz);
    sum_axpy(- beta * omega, v0, p0, nx, ny, nz);
    pixel_data pv(boost::extents[n_angles][n_h][n_v]);
    device->forward_project(pv, p0, origin, voxel_size,
			    (int)sz[0], (int)sz[1], (int)sz[2]);
    init_data(v0, nx, ny, nz);
    device->backward_project(v0, pv, origin, voxel_size,
			     (int)sz[0], (int)sz[1], (int)sz[2]);
    alpha = gamma / norm_voxels(r0, v0, nx, ny, nz);
    update_progress(2 * iter + 2);
    voxel_data s(boost::extents[sz[0]][sz[1]][sz[2]],
		 boost::c_storage_order());
    // Todo - combine these
    copy(r, s, nx, ny, nz);
    sum_axpy(- alpha, v0, s, nx, ny, nz);
    init_data(pv, n_angles, n_h, n_v);
    device->forward_project(pv, s, origin, voxel_size,
			    (int)sz[0], (int)sz[1], (int)sz[2]);
    voxel_data t(boost::extents[sz[0]][sz[1]][sz[2]],
		 boost::c_storage_order());
    device->backward_project(t, pv, origin, voxel_size,
			     (int)sz[0], (int)sz[1], (int)sz[2]);
    omega = norm_voxels(t, s, nx, ny, nz) / norm_voxels(t, nx, ny, nz);
    sum_axpy(alpha, p0, voxels, nx, ny, nz);
    sum_axpy(omega, s, voxels, nx, ny, nz);
    // Todo - combine
    copy(s, r, nx, ny, nz);
    sum_axpy(- omega, t, r, nx, ny, nz);
    gamma0  = gamma;
    update_progress(2 * iter + 3);
    iter_time.accumulate();
    iter_time.output("Iteration ");
  }
  return true;
}
