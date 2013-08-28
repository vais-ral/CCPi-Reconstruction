
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
  const real ninety = M_PI / 2.0;
  const real one_eighty = M_PI;
  const real two_seventy = 3.0 * M_PI / 2.0;
  const real three_sixty = 2.0 * M_PI;
  //real xmin = vox_origin[0];
  //real xmax = vox_origin[0] + real(nx) * vox_size[0];  
  //real ymin = vox_origin[1];
  //real ymax = vox_origin[1] + real(ny) * vox_size[1];
  //real extent = std::max(std::max(std::abs(xmin), std::abs(xmax)),
  //			 std::max(std::abs(ymin), std::abs(ymax)));
  // diagonal of slice is sqrt(2) (max - min), so use 2 * extent should be ok
  //real p = 2.0 * extent;
  std::vector<real> cangles(n_angles);
  for (int i = 0; i < n_angles; i++)
    cangles[i] = std::cos(angles[i]);
  std::vector<real> sangles(n_angles);
  for (int i = 0; i < n_angles; i++)
    sangles[i] = std::sin(angles[i]);
  std::vector<real> x_offsets(nx + 1);
  for (int i = 0; i <= nx; i++)
    x_offsets[i] = real(i) * vox_size[0];
  std::vector<real> x_coords(nx + 1);
  for (int i = 0; i <= nx; i++)
    x_coords[i] = vox_origin[0] + x_offsets[i];
  std::vector<real> y_coords(ny + 1);
  for (int i = 0; i <= ny; i++)
    y_coords[i] = vox_origin[1] + real(i) * vox_size[1];
  std::vector<real> z_coords(nz + 1);
  for (int i = 0; i <= nz; i++)
    z_coords[i] = vox_origin[2] + real(i) * vox_size[2];
  std::vector<real> line_lengths(n_angles);
  for (int i = 0; i < n_angles; i++) {
    if (std::abs(sangles[i]) < 1e-6)
      line_lengths[i] = 1e8;
    else
      line_lengths[i] = vox_size[1] / sangles[i];
  }
  real pixel_step = h_pixels[1] - h_pixels[0];
  // Todo - build an array of rotated pixel px/py?
  std::vector<real> vox_lengths(n_angles);
  for (int i = 0; i < n_angles; i++)
    vox_lengths[i] = vox_size[0] / cangles[i];
  real inv_vox_x = 1.0 / vox_size[0];
  //real inv_vox_y = 1.0 / vox_size[1];
  //real inv_vox_z = 1.0 / vox_size[2];
  // put end value to stop while loop without over-run
  std::vector<real> v_coords(nv_pixels + 1);
  for (int i = 0; i < nv_pixels; i++)
    v_coords[i] = v_pixels[i];
  v_coords[nv_pixels] = z_coords[nz] + 1.0;
  int v = 0;
  while (v_pixels[v] < z_coords[0])
    v++;
  for (int k = 0; k < nz; k++) {
    while (v_coords[v] < z_coords[k + 1]) {
      // Todo - what is best order of a/j loops and which to parallelise?
      for (int a = 0; a < n_angles; a++) {
	real phi = angles[a];
	real cphi = cangles[a];
	real sphi = sangles[a];
	real LL = line_lengths[a];
	// something about > 45 and L > LL to worry about below?
	// especially if i1 + 1 == i2 and 2 calcs not part of loop? Todo
	real L = vox_lengths[a];
	if (std::abs(sphi) < 1e-5) {
	  // 0.0 or 180.0
	  for (int j = 0; j < ny; j++) {
	    long vox_offset = long(k) * nx * ny + long(j) * nx;
	    real y1 = y_coords[j];
	    real y2 = y_coords[j + 1];
	    for (int h = 0; h < nh_pixels; h++) {
	      long pixel_index = long(a) * nv_pixels * nh_pixels
		+ long(v) * nh_pixels + long(h);
	      real py = h_pixels[h];
	      // only intercepts if within boundary and then intecepts all
	      if (y1 <= py and py < y2) {
		for (int i = 0; i < nx; i++)
		  voxels[vox_offset + i] += pixels[pixel_index] * L;
	      }
	    }
	  }
	} else if (phi > 0.0 and phi <= ninety) {
	  real s_scale = 1.0 / sphi;
	  for (int j = 0; j < ny; j++) {
	    long vox_offset = long(k) * nx * ny + long(j) * nx;
	    real y1 = y_coords[j];
	    real y2 = y_coords[j + 1];
	    int h0 = int((y2 * cphi - h_pixels[0]
			  - x_coords[0] * sphi) / pixel_step);
	    int hN = int((y1 * cphi - h_pixels[0]
			  - x_coords[nx] * sphi) / pixel_step);
	    // precalc y0/yn values?
	    if (hN < 0)
	      hN = 0;
	    if (h0 >= nh_pixels)
	      h0 = nh_pixels;
	    //std::cout << "range " << hN << ' ' << h0 << '\n';
	    for (int h = hN; h < h0; h++) {
	      long pixel_index = long(a) * nv_pixels * nh_pixels
		+ long(v) * nh_pixels + long(h);
	      real py = h_pixels[h];
	      real r1 = (y1 - py * cphi) * s_scale;
	      // Todo - surely this is fixed for each pixel step?
	      //real rstep = (vox_size[1] - py * cphi) * s_scale;
	      //real r2 = r1 + rstep;
	      real r2 = (y2 - py * cphi) * s_scale;
	      //real rstep = r2 - r1;
	      real x1 = - py * sphi + r1 * cphi;
	      real x2 = - py * sphi + r2 * cphi;
	      // range of i for strip
	      // Todo - check these truncate the way I expect
	      int i1 = (x1 - vox_origin[0]) * inv_vox_x;
	      int i2 = (x2 - vox_origin[0]) * inv_vox_x;
	      if (i1 <= i2 and i2 >= 0 and i1 < nx) {
		if (i1 == i2)
		  voxels[vox_offset + i1] += pixels[pixel_index] * LL;
		else if (i2 == 0) {
		  voxels[vox_offset + i2] += pixels[pixel_index]
		    * (x2 - x_coords[0]) * inv_vox_x * L;
		} else if (i1 == nx - 1) {
		  voxels[vox_offset + i1] += pixels[pixel_index]
		    * (x_coords[nx] - x1) * inv_vox_x * L;
		} else {
		  if (i1 < 0)
		    i1 = 0;
		  if (i2 >= nx)
		    i2 = nx - 1;
		  //real Lextra = LL - real(i2 - i1 - 1) * L;
		  // Todo - simplify diff from x1 as both use origin
		  real xs = x_coords[i1 + 1];
		  // ratio of base used to length if full  base used
		  voxels[vox_offset + i1] += pixels[pixel_index]
		    * (xs - x1) * inv_vox_x * L;
		  for (int i = i1 + 1; i < i2; i++)
		    voxels[vox_offset + i] += pixels[pixel_index] * L;
		  real xe = x_coords[i2];
		  // ratio of base used to length if full  base used
		  voxels[vox_offset + i2] += pixels[pixel_index]
		    * (x2 - xe) * inv_vox_x * L;
		}
	      }
	    }
	  }
	} else if (phi > ninety and phi < one_eighty) {
	  real s_scale = 1.0 / sphi;
	  for (int j = 0; j < ny; j++) {
	    long vox_offset = long(k) * nx * ny + long(j) * nx;
	    real y1 = y_coords[j];
	    real y2 = y_coords[j + 1];
	    // precalc y0/yn values?
	    for (int h = 0; h < nh_pixels; h++) {
	      long pixel_index = long(a) * nv_pixels * nh_pixels
		+ long(v) * nh_pixels + long(h);
	      real py = h_pixels[h];
	      real r1 = (y1 - py * cphi) * s_scale;
	      // Todo - surely this is fixed for each pixel step?
	      //real rstep = (vox_size[1] - py * cphi) * s_scale;
	      //real r2 = r1 + rstep;
	      real r2 = (y2 - py * cphi) * s_scale;
	      //real rstep = r2 - r1;
	      real x1 = - py * sphi + r1 * cphi;
	      real x2 = - py * sphi + r2 * cphi;
	      // range of i for strip
	      // Todo - check these truncate the way I expect
	      int i1 = (x1 - vox_origin[0]) * inv_vox_x;
	      int i2 = (x2 - vox_origin[0]) * inv_vox_x;
	      // Todo need to block these? i1 = i2 or i1 +1 = i2 only??
	      if (i2 <= i1 and i1 >= 0 and i2 < nx) {
		if (i1 == i2)
		  voxels[vox_offset + i1] += pixels[pixel_index] * LL;
		else if (i1 == 0)
		  voxels[vox_offset + i1] += pixels[pixel_index]
		    * (x1 - x_coords[0]) * inv_vox_x * L;
		else if (i2 == nx - 1)
		  voxels[vox_offset + i2] += pixels[pixel_index]
		    * (x_coords[nx] - x2) * inv_vox_x * L;
		else {
		  if (i2 < 0)
		    i2 = 0;
		  if (i1 >= nx)
		    i1 = nx - 1;
		  //real Lextra = LL - real(i1 - i2 - 1) * L;
		  // Todo - simplify diff from x1 as both use origin
		  real xs = x_coords[i2 + 1];
		  // ratio of base used to length if full  base used
		  voxels[vox_offset + i2] += pixels[pixel_index]
		    * (xs - x2) * inv_vox_x * L;
		  for (int i = i2 + 1; i < i1; i++)
		    voxels[vox_offset + i] += pixels[pixel_index] * L;
		  real xe = x_coords[i1];
		  // ratio of base used to length if full  base used
		  voxels[vox_offset + i1] += pixels[pixel_index]
		    * (x1 - xe) * inv_vox_x * L;
		}
	      }
	    }
	  }
	} else if (phi > one_eighty and phi < two_seventy) {
	  real s_scale = 1.0 / sphi;
	  for (int j = 0; j < ny; j++) {
	    long vox_offset = long(k) * nx * ny + long(j) * nx;
	    real y1 = y_coords[j];
	    real y2 = y_coords[j + 1];
	    // precalc y0/yn values?
	    for (int h = 0; h < nh_pixels; h++) {
	      long pixel_index = long(a) * nv_pixels * nh_pixels
		+ long(v) * nh_pixels + long(h);
	      real py = h_pixels[h];
	      real r1 = (y1 - py * cphi) * s_scale;
	      // Todo - surely this is fixed for each pixel step?
	      //real rstep = (vox_size[1] - py * cphi) * s_scale;
	      //real r2 = r1 + rstep;
	      real r2 = (y2 - py * cphi) * s_scale;
	      //real rstep = r2 - r1;
	      real x1 = - py * sphi + r1 * cphi;
	      real x2 = - py * sphi + r2 * cphi;
	      // range of i for strip
	      // Todo - check these truncate the way I expect
	      int i1 = (x1 - vox_origin[0]) * inv_vox_x;
	      int i2 = (x2 - vox_origin[0]) * inv_vox_x;
	      // Todo need to block these? i1 = i2 or i1 +1 = i2 only??
	      if (i2 <= i1 and i1 >= 0 and i2 < nx) {
		if (i1 == i2)
		  voxels[vox_offset + i1] += pixels[pixel_index] * LL;
		else if (i1 == 0)
		  voxels[vox_offset + i1] += pixels[pixel_index]
		    * (x1 - x_coords[0]) * inv_vox_x * L;
		else if (i2 == nx - 1)
		  voxels[vox_offset + i2] += pixels[pixel_index]
		    * (x_coords[nx] - x2) * inv_vox_x * L;
		else {
		  if (i2 < 0)
		    i2 = 0;
		  if (i1 >= nx)
		    i1 = nx - 1;
		  //real Lextra = LL - real(i1 - i2 - 1) * L;
		  // Todo - simplify diff from x1 as both use origin
		  real xs = x_coords[i2 + 1];
		  // ratio of base used to length if full  base used
		  voxels[vox_offset + i2] += pixels[pixel_index]
		    * (xs - x2) * inv_vox_x * L;
		  for (int i = i2 + 1; i < i1; i++)
		    voxels[vox_offset + i] += pixels[pixel_index] * L;
		  real xe = x_coords[i1];
		  // ratio of base used to length if full  base used
		  voxels[vox_offset + i1] += pixels[pixel_index]
		    * (x1 - xe) * inv_vox_x * L;
		}
	      }
	    }
	  }
	} else if (phi >= two_seventy and phi < three_sixty) {
	  real s_scale = 1.0 / sphi;
	  for (int j = 0; j < ny; j++) {
	    long vox_offset = long(k) * nx * ny + long(j) * nx;
	    real y1 = y_coords[j];
	    real y2 = y_coords[j + 1];
	    // precalc y0/yn values?
	    for (int h = 0; h < nh_pixels; h++) {
	      long pixel_index = long(a) * nv_pixels * nh_pixels
		+ long(v) * nh_pixels + long(h);
	      real py = h_pixels[h];
	      real r1 = (y1 - py * cphi) * s_scale;
	      // Todo - surely this is fixed for each pixel step?
	      //real rstep = (vox_size[1] - py * cphi) * s_scale;
	      //real r2 = r1 + rstep;
	      real r2 = (y2 - py * cphi) * s_scale;
	      //real rstep = r2 - r1;
	      real x1 = - py * sphi + r1 * cphi;
	      real x2 = - py * sphi + r2 * cphi;
	      // range of i for strip
	      // Todo - check these truncate the way I expect
	      int i1 = (x1 - vox_origin[0]) * inv_vox_x;
	      int i2 = (x2 - vox_origin[0]) * inv_vox_x;
	      if (i1 <= i2 and i2 >= 0 and i1 < nx) {
		if (i1 == i2)
		  voxels[vox_offset + i1] += pixels[pixel_index] * LL;
		else if (i2 == 0) {
		  voxels[vox_offset + i2] += pixels[pixel_index]
		    * (x2 - x_coords[0]) * inv_vox_x * L;
		} else if (i1 == nx - 1) {
		  voxels[vox_offset + i1] += pixels[pixel_index]
		    * (x_coords[nx] - x1) * inv_vox_x * L;
		} else {
		  if (i1 < 0)
		    i1 = 0;
		  if (i2 >= nx)
		    i2 = nx - 1;
		  //real Lextra = LL - real(i2 - i1 - 1) * L;
		  // Todo - simplify diff from x1 as both use origin
		  real xs = x_coords[i1 + 1];
		  // ratio of base used to length if full  base used
		  voxels[vox_offset + i1] += pixels[pixel_index]
		    * (xs - x1) * inv_vox_x * L;
		  for (int i = i1 + 1; i < i2; i++)
		    voxels[vox_offset + i] += pixels[pixel_index] * L;
		  real xe = x_coords[i2];
		  // ratio of base used to length if full  base used
		  voxels[vox_offset + i2] += pixels[pixel_index]
		    * (x2 - xe) * inv_vox_x * L;
		}
	      }
	    }
	  }
	} else
	  std::cerr << "Angle out of range\n";
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
    /*
    instrument::backward_project(get_h_pixels(), get_v_pixels(), get_phi(),
				 get_pixel_data(), voxels,
				 get_num_angles(), get_num_h_pixels(),
				 get_num_v_pixels(), origin, width, nx, ny, nz);
  */
    my_back_project(get_h_pixels(), get_v_pixels(), get_phi(),
				 get_pixel_data(), voxels,
				 get_num_angles(), get_num_h_pixels(),
				 get_num_v_pixels(), origin, width, nx, ny, nz);
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
