
#include <cmath>
#include <iostream>
#include <omp.h>
#ifdef MATLAB_MEX_FILE
#  include "mex_types.hpp"
#else
#  include "base_types.hpp"
#endif // mex
#include "instruments.hpp"
#include "project_line.hpp"
#include "parallel_b.hpp"
#include "parallel_f.hpp"
#include "cone_b.hpp"
#include "cone_f.hpp"
#include "timer.hpp"

#ifdef MKL_ILP64
#  define MKL_INT long
#  include "mkl_spblas.h"
#endif // MKL_ILP64

#ifndef USE_TIMER
#  define USE_TIMER false
#endif // USE_TIMER

/*
   Todo - should improve structure of pixel_data and not require the rest
   of the code to know the order, although this may be hard to achieve with the
   current Matlab structures.
*/

long CCPi::instrument::get_data_size() const
{
  return long(n_angles) * long(n_vertical_pixels) * long(n_horizontal_pixels);
}

pixel_type *const CCPi::instrument::get_pixel_data() const
{
  return pixel_data;
}

pixel_type *CCPi::instrument::create_pixel_data()
{
  long n_rays = long(n_angles) * long(n_vertical_pixels)
    * long(n_horizontal_pixels);
  pixel_data = new pixel_type[n_rays];
  return pixel_data;
}

void CCPi::instrument::set_pixel_data(pixel_type *p, const long n)
{
  long n_rays = long(n_angles) * long(n_vertical_pixels)
    * long(n_horizontal_pixels);
  if (n != n_rays)
    std::cerr << "Size mismatch setting pixel data\n";
  pixel_data = p;
}

void CCPi::instrument::set_phi(real *p, const int n)
{
  phi = p;
  n_angles = n;
}

void CCPi::cone_beam::set_params(const real sx, const real sy, const real sz,
				 const real dx, real dy[], real dz[],
				 real ang[], const int ny, const int nz,
				 const int nang)
{
  set_source(sx, sy, sz);
  set_detector(dx);
  set_h_pixels(dy, ny);
  set_v_pixels(dz, nz);
  set_phi(ang, nang);
}

void CCPi::cone_beam::forward_project(pixel_type *pixels,
				      voxel_type *const voxels,
				      const real origin[3],
				      const real width[3], const int nx,
				      const int ny, const int nz) const
{
  timer fptime(USE_TIMER);
  instrument::forward_project(source_x, source_y, source_z, detector_x,
			      get_h_pixels(), get_v_pixels(), get_phi(),
			      pixels, voxels, get_num_angles(),
			      get_num_h_pixels(), get_num_v_pixels(), origin,
			      width, nx, ny, nz);
  fptime.accumulate();
  fptime.output(" forward projection");
}

void CCPi::cone_beam::backward_project(pixel_type *pixels,
				       voxel_type *const voxels,
				       const real origin[3],
				       const real width[3], const int nx,
				       const int ny, const int nz) const
{
  timer bptime(USE_TIMER);
  instrument::backward_project(source_x, source_y, source_z, detector_x,
			       get_h_pixels(), get_v_pixels(), get_phi(),
			       pixels, voxels, get_num_angles(),
			       get_num_h_pixels(), get_num_v_pixels(), origin,
			       width, nx, ny, nz);
  bptime.accumulate();
  bptime.output("backward projection");
}

void CCPi::cone_beam::backward_project(voxel_type *const voxels,
				       const real origin[3],
				       const real width[3], const int nx,
				       const int ny, const int nz) const
{
  timer bptime(USE_TIMER);
  instrument::backward_project(source_x, source_y, source_z, detector_x,
			       get_h_pixels(), get_v_pixels(), get_phi(),
			       get_pixel_data(), voxels,
			       get_num_angles(), get_num_h_pixels(),
			       get_num_v_pixels(), origin, width, nx, ny, nz);
  bptime.accumulate();
  bptime.output("backward projection");
}

void CCPi::parallel_beam::forward_project(pixel_type *pixels,
					  voxel_type *const voxels,
					  const real origin[3],
					  const real width[3], const int nx,
					  const int ny, const int nz) const
{
  timer fptime(USE_TIMER);
  if (has_projection_matrix)
    forward_project_matrix(get_v_pixels(), pixels, voxels, get_num_angles(),
			   get_num_h_pixels(), get_num_v_pixels(), origin,
			   width, nx, ny, nz);
  else
    instrument::forward_project(get_h_pixels(), get_v_pixels(), get_phi(),
				pixels, voxels, get_num_angles(),
				get_num_h_pixels(), get_num_v_pixels(), origin,
				width, nx, ny, nz);
  fptime.accumulate();
  fptime.output(" forward projection");
}

extern void test_voxel(const int i, const int j, const int k,
		       const int a, const int v, const int h, const double d,
		       const double l, const double voxel_origin[3],
		       const double voxel_size[3], const double v_pos,
		       const double h_pos, const double cphi,
		       const double sphi, const int nx, const int ny,
		       const int nz, const char message[]);

// Todo - simdise across a group of j or h?
void CCPi::parallel_beam::my_back_project(const real h_pixels[],
					  const real v_pixels[],
					  const real angles[],
					  pixel_type pixels[],
					  voxel_type *const voxels,
					  const int n_angles,
					  const int nh_pixels,
					  const int nv_pixels,
					  const real vox_origin[3],
					  const real vox_size[3],
					  const int nx,
					  const int ny,
					  const int nz)
{
  //const real tol = 1e-13;
  std::vector<real> cangles(n_angles);
  for (int i = 0; i < n_angles; i++)
    cangles[i] = std::cos(angles[i]);
  std::vector<real> sangles(n_angles);
  for (int i = 0; i < n_angles; i++)
    sangles[i] = std::sin(angles[i]);
  // centred positions
  std::vector<real> x0_coords(nx + 1);
  for (int i = 0; i <= nx; i++)
    x0_coords[i] = vox_origin[0] + real(i) * vox_size[0];
  std::vector<real> y0_coords(ny + 1);
  for (int i = 0; i <= ny; i++)
    y0_coords[i] = vox_origin[1] + real(i) * vox_size[1];
  std::vector<real> z_coords(nz + 1);
  for (int i = 0; i <= nz; i++)
    z_coords[i] = vox_origin[2] + real(i) * vox_size[2];
  real inv_pixel_step = 1.0 / (h_pixels[1] - h_pixels[0]);
  // Todo - pre-rotated array of voxels? need for each angle, x, y though
  // x' = x cos() + y sin()
  // y' = -x sin() + y cos() - not used in current form
  std::vector<real> vox_lengths(n_angles);
  for (int i = 0; i < n_angles; i++) {
    if (std::abs(cangles[i]) > std::abs(sangles[i]))
      vox_lengths[i] = vox_size[0] / std::abs(cangles[i]);
    else
      vox_lengths[i] = vox_size[0] / std::abs(sangles[i]);
  }
  // Todo can precalc 2D i0/inh1, yc0, xh_coords here
  // its probably p that is the tricky one though
  // put end value to stop while loop without over-run
  std::vector<real> v_coords(nv_pixels + 1);
  for (int i = 0; i < nv_pixels; i++)
    v_coords[i] = v_pixels[i];
  v_coords[nv_pixels] = z_coords[nz] + 1.0;
  /*std::vector<real> lengths(nh_pixels);*/
  std::vector<real> xx_coords(nx + 1);
  std::vector<real> yy_coords(ny + 1);
  int v = 0;
  while (v_pixels[v] < z_coords[0])
    v++;
  for (int k = 0; k < nz; k++) {
    long vox_z = k * long(ny) * long(nx);
    // make this an inner loop? Todo
    while (v_coords[v] < z_coords[k + 1]) {
      long pix_v = v * long(nh_pixels);
      // Todo - what is best order of a/j loops and which to parallelise?
      for (int a = 0; a < n_angles; a++) {
	long pix_av = a * long(nv_pixels) * long(nh_pixels) + pix_v;
	real cphi = cangles[a];
	real sphi = sangles[a];
	real L = vox_lengths[a];
	//for (int i = 0; i < nx; i++)
	//xh_coords[i] = - x_coords[i] * sphi * inv_pixel_step;
	for (int i = 0; i <= nx; i++)
	  xx_coords[i] = - x0_coords[i] * sphi;
	for (int i = 0; i <= ny; i++)
	  yy_coords[i] = y0_coords[i] * cphi;
	real y_00 = xx_coords[0] + yy_coords[0];
	real y_01 = xx_coords[0] + yy_coords[1];
	real y_10 = xx_coords[1] + yy_coords[0];
	real y_11 = xx_coords[1] + yy_coords[1];
	real y_top = 0.0;
	real y_l1;
	real y_l2;
	real y_bot = 0.0;
	if (sphi >= 0.0) {
	  // <= 180, scanning from 0 to nx decreases p
	  // orientation of voxels is the same for all so can get rotated shape
	  // from any single one - use 0.0
	  int jm;
	  if (cphi >= 0.0) { // 0-90
	    jm = 1;
	    if (cphi >= sphi) { // 0-45
	      y_l1 = y_01 - y_11;
	      y_l2 = y_01 - y_00;
	      y_bot = y_01 - y_10;
	    } else {
	      y_l1 = y_01 - y_00;
	      y_l2 = y_01 - y_11;
	      y_bot = y_01 - y_10;
	    }
	  } else {
	    jm = 0;
	    if (std::abs(cphi) > sphi) { // 135-180
	      y_l1 = y_00 - y_10;
	      y_l2 = y_00 - y_01;
	      y_bot = y_00 - y_11;
	    } else {
	      y_l1 = y_00 - y_01;
	      y_l2 = y_00 - y_10;
	      y_bot = y_00 - y_11;
	    }
	  }
	  //const int npixels = int(y_bot * inv_pixel_step) + 2;
	  // Its always scaled by L
	  real inv_y_l1 = L / y_l1;
	  for (int j = 0; j < ny; j++) {
	    real yj = yy_coords[j + jm];
	    long vox_zy = vox_z + j * long(nx);
	    int p = int((xx_coords[0] + yj - h_pixels[0]) * inv_pixel_step);
	    if (p >= nh_pixels)
	      p = nh_pixels - 1;
	    else if (p < 0)
	      continue;
	    for (int i = 0; i < nx; i++) {
	      // Todo ? ytop decreases by y_l1 for each i step
	      y_top = xx_coords[i] + yj;
	      while (h_pixels[p] > y_top) {
		p--;
		if (p < 0)
		  break;
	      }
	      if (p < 0)
		break;
	      real yb = y_top - y_bot;
	      real y1 = y_top - y_l1;
	      real y2 = y_top - y_l2;
	      /*int cnt = 0;*/
	      for (int q = p; q > -1; q--) {
		if (h_pixels[q] < yb)
		  break;
		else if (h_pixels[q] < y2) {
		  real ratio = (h_pixels[q] - yb) * inv_y_l1;
		  /*lengths[cnt] = ratio;*/
		  voxels[vox_zy + i] += pixels[pix_av + q] * ratio;
		} else if (h_pixels[q] < y1) {
		  /*lengths[cnt] = L;*/
		  voxels[vox_zy + i] += pixels[pix_av + q] * L;
		} else {
		  real ratio = (y_top - h_pixels[q]) * inv_y_l1;
		  /*lengths[cnt] = ratio;*/
		  voxels[vox_zy + i] += pixels[pix_av + q] * ratio;
		}
		/*cnt++;*/
	      }
	      /*
	      if (k == 0) {
		if (p < nh_pixels - 1)
		  test_voxel(i, j, k, a, v, p+1, 0.0, 0.0,
			     vox_origin, vox_size, v_coords[v], h_pixels[p+1],
			     cphi, sphi, nx, ny, nz, "p+");
		for (int q = 0; q < cnt; q++) {
		  test_voxel(i, j, k, a, v, p-q, 1.0, lengths[q],
			     vox_origin, vox_size, v_coords[v], h_pixels[p-q],
			     cphi, sphi, nx, ny, nz, "pq");
		}
		if (p-cnt >= 0) {
		  test_voxel(i, j, k, a, v, p-cnt, 0.0, 0.0,
			     vox_origin, vox_size, v_coords[v], h_pixels[p-cnt],
			     cphi, sphi, nx, ny, nz, "p-");
		}
	      }
	      */
	    }
	  }
	} else {
	  // 0-nx increases p
	  int jm;
	  if (cphi < 0.0) { // 180-270
	    jm = 1;
	    if (cphi < sphi) { // 180-225
	      y_l1 = y_11 - y_01;
	      y_l2 = y_00 - y_01;
	      y_top = y_10 - y_01;
	    } else {
	      y_l1 = y_00 - y_01;
	      y_l2 = y_11 - y_01;
	      y_top = y_10 - y_01;
	    }
	  } else { //270-360
	    jm = 0;
	    if (cphi > std::abs(sphi)) { // 315-360
	      y_l1 = y_10 - y_00;
	      y_l2 = y_01 - y_00;
	      y_top = y_11 - y_00;
	    } else {
	      y_l1 = y_01 - y_00;
	      y_l2 = y_10 - y_00;
	      y_top = y_11 - y_00;
	    }
	  }
	  real inv_y_l1 = L / y_l1;
	  for (int j = 0; j < ny; j++) {
	    real yj = yy_coords[j + jm];
	    long vox_zy = vox_z + j * long(nx);
	    int p = int((xx_coords[0] + yj - h_pixels[0]) * inv_pixel_step);
	    if (p < 0)
	      p = 0;
	    else if (p >= nh_pixels)
	      continue;
	    for (int i = 0; i < nx; i++) {
	      // Todo ? ybot increases by y_l1 for each i step
	      y_bot = xx_coords[i] + yj;
	      while (h_pixels[p] < y_bot) {
		p++;
		if (p >= nh_pixels)
		  break;
	      }
	      if (p >= nh_pixels)
		break;
	      real yt = y_bot + y_top;
	      real y1 = y_bot + y_l1;
	      real y2 = y_bot + y_l2;
	      /*int cnt = 0;*/
	      for (int q = p; q < nh_pixels; q++) {
		if (h_pixels[q] >= yt)
		  break;
		else if (h_pixels[q] < y1) {
		  real ratio = (h_pixels[q] - y_bot) * inv_y_l1;
		  /*lengths[cnt] = ratio;*/
		  voxels[vox_zy + i] += pixels[pix_av + q] * ratio;
		} else if (h_pixels[q] < y2) {
		  /*lengths[cnt] = L;*/
		  voxels[vox_zy + i] += pixels[pix_av + q] * L;
		} else {
		  real ratio = (yt - h_pixels[q]) * inv_y_l1;
		  /*lengths[cnt] = ratio;*/
		  voxels[vox_zy + i] += pixels[pix_av + q] * ratio;
		}
		/*cnt++;*/
	      }
	      /*
	      if (k == 0) {
		if (i == 0 and j == 0 and a == 46 and v == 0) {
		  std::cerr << p << ' ' << cnt << ' ' << y_top << ' '
			    << y_l1 << ' ' << y_l2 << ' ' << y_bot << ' '
			    << h_pixels[p] << '\n';
		}
		if (p > 0)
		  test_voxel(i, j, k, a, v, p-1, 0.0, 0.0,
			     vox_origin, vox_size, v_coords[v], h_pixels[p-1],
			     cphi, sphi, nx, ny, nz, "p-");
		for (int q = 0; q < cnt; q++) {
		  test_voxel(i, j, k, a, v, p+q, 1.0, lengths[q],
			     vox_origin, vox_size, v_coords[v], h_pixels[p+q],
			     cphi, sphi, nx, ny, nz, "pq");
		}
		if (p+cnt < nh_pixels) {
		  test_voxel(i, j, k, a, v, p+cnt, 0.0, 0.0,
			     vox_origin, vox_size, v_coords[v], h_pixels[p+cnt],
			     cphi, sphi, nx, ny, nz, "p+");
		}
	      }
	      */
	    }
	  }
	}
      }
      v++;
    }
  }
}

void CCPi::parallel_beam::backward_project(pixel_type *pixels,
					   voxel_type *const voxels,
					   const real origin[3],
					   const real width[3], const int nx,
					   const int ny, const int nz) const
{
  timer bptime(USE_TIMER);
  if (has_projection_matrix)
    backward_project_matrix(get_v_pixels(), pixels, voxels, get_num_angles(),
			    get_num_h_pixels(), get_num_v_pixels(), origin,
			    width, nx, ny, nz);
  else
    /*
    instrument::backward_project(get_h_pixels(), get_v_pixels(), get_phi(),
				 pixels, voxels, get_num_angles(),
				 get_num_h_pixels(), get_num_v_pixels(), origin,
				 width, nx, ny, nz);
    */
    my_back_project(get_h_pixels(), get_v_pixels(), get_phi(),
				 pixels, voxels, get_num_angles(),
				 get_num_h_pixels(), get_num_v_pixels(), origin,
				 width, nx, ny, nz);
  bptime.accumulate();
  bptime.output("backward projection");
}

void CCPi::parallel_beam::backward_project(voxel_type *const voxels,
					   const real origin[3],
					   const real width[3], const int nx,
					   const int ny, const int nz) const
{
  timer bptime(USE_TIMER);
  if (has_projection_matrix)
    backward_project_matrix(get_v_pixels(), get_pixel_data(), voxels,
			    get_num_angles(), get_num_h_pixels(),
			    get_num_v_pixels(), origin, width, nx, ny, nz);
  else
    /**/
    instrument::backward_project(get_h_pixels(), get_v_pixels(), get_phi(),
				 get_pixel_data(), voxels,
				 get_num_angles(), get_num_h_pixels(),
				 get_num_v_pixels(), origin, width, nx, ny, nz);
  /**/ /*
    my_back_project(get_h_pixels(), get_v_pixels(), get_phi(),
				 get_pixel_data(), voxels,
				 get_num_angles(), get_num_h_pixels(),
				 get_num_v_pixels(), origin, width, nx, ny, nz);
       */
  bptime.accumulate();
  bptime.output("backward projection");
}

void CCPi::cone_beam::setup_projection_matrix(const real origin[3],
					      const real width[3],
					      const int nx, const int ny,
					      const int nz)
{
  std::cerr << "Projection matrix not implemented\n";
}

void CCPi::parallel_beam::setup_projection_matrix(const real origin[3],
						  const real width[3],
						  const int nx, const int ny,
						  const int nz)
{
  // Should really check that the sizes are the same.
  if (has_projection_matrix)
    return;
  timer ptime(USE_TIMER);
  // assumes 2D slices which makes the storage practical.
  setup_2D_matrix(get_h_pixels(), get_phi(), get_num_angles(),
		  get_num_v_pixels(), get_num_h_pixels(), origin, width,
		  nx, ny, nz);
  has_projection_matrix = true;
  ptime.accumulate();
  ptime.output("projection map");
}

void CCPi::parallel_beam::forward_project_matrix(const real det_z[],
						 pixel_type ray_data[],
						 voxel_type *const vol_data,
						 const int n_angles,
						 const int n_rays_y,
						 const int n_rays_z,
						 const real grid_offset[3],
						 const real voxel_size[3],
						 const int nx_voxels,
						 const int ny_voxels,
						 const int nz_voxels) const
{
#pragma omp parallel for shared(det_z, ray_data) schedule(dynamic)
  for (long curr_ray_z = 0; curr_ray_z < n_rays_z; curr_ray_z++) {
    //real z = det_z[curr_ray_z];
    long k = (long)std::floor((det_z[curr_ray_z]
			       - grid_offset[2]) / voxel_size[2]);
    if (k < 0 or k >= nz_voxels)
      continue;
    long k_offset = k * nx_voxels * ny_voxels;
    long z_offset = curr_ray_z * n_rays_y;

    for (long i = 0; i < matrix_size; i++)
      ray_data[forward_rows[i] + z_offset] += forward_matrix[i]
	* vol_data[k_offset + forward_cols[i]];
    /*
    char desc[6];
    desc[0] = 'G';
    desc[3] = 'C';
    long n_rays = long(n_angles) * long(n_rays_z) * long(n_rays_y);
    long n_vox = long(nx_voxels) * long(ny_voxels) * long(nz_voxels);
    real alpha = 1.0;
    real beta = 1.0;
    long sz = matrix_size;
    mkl_dcoomv("N", &n_rays, &n_vox, &alpha, desc, forward_matrix,
	       forward_rows, forward_cols, &sz, &vol_data[k_offset],
	       &beta, &ray_data[z_offset]);
    */
  }
}

void CCPi::parallel_beam::backward_project_matrix(const real det_z[],
						  pixel_type ray_data[],
						  voxel_type *const vol_data,
						  const int n_angles,
						  const int n_rays_y,
						  const int n_rays_z,
						  const real grid_offset[3],
						  const real voxel_size[3],
						  const int nx_voxels,
						  const int ny_voxels,
						  const int nz_voxels) const
{
#pragma omp parallel for shared(det_z, ray_data) schedule(dynamic)
  for (long k = 0; k < nz_voxels; k++) {
    long k_offset = k * nx_voxels * ny_voxels;
    long z_min = -1;
    long z_max = -1;
    // Todo - calculate range rather than looping
    for (long curr_ray_z = 0; curr_ray_z < n_rays_z; curr_ray_z++) {
      long kz = (long)std::floor((det_z[curr_ray_z]
				  - grid_offset[2]) / voxel_size[2]);
      if (kz == k) {
	if (z_min == -1)
	  z_min = curr_ray_z;
	z_max = curr_ray_z;
      } else if (kz > k)
	break;
    }

    long size = nx_voxels * ny_voxels;
    for (long curr_ray_z = z_min; curr_ray_z <= z_max; curr_ray_z++) {
      long z_offset = curr_ray_z * n_rays_y;
      for (long i = 0; i < size; i++) {
	for (long j = backward_rowb[i]; j < backward_rowe[i]; j++) {
	  vol_data[k_offset + i] += backward_matrix[j] *
	    ray_data[backward_cols[j] + z_offset];
	}
      }
      /*
      char desc[6];
      desc[0] = 'G';
      desc[3] = 'C';
      long n_rays = long(n_angles) * long(n_rays_z) * long(n_rays_y);
      real alpha = 1.0;
      real beta = 1.0;
      mkl_dcsrmv("N", &size, &n_rays, &alpha, desc, backward_matrix,
		 backward_cols, backward_rowb, backward_rowe,
		 &ray_data[z_offset], &beta, &vol_data[k_offset]);
      */
    }
  }
}
