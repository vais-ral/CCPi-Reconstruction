
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

sl_int CCPi::instrument::get_data_size() const
{
  return sl_int(n_angles) * sl_int(n_vertical_pixels)
    * sl_int(n_horizontal_pixels);
}

pixel_type *const CCPi::instrument::get_pixel_data() const
{
  return pixel_data;
}

pixel_type *CCPi::instrument::create_pixel_data()
{
  sl_int n_rays = sl_int(n_angles) * sl_int(n_vertical_pixels)
    * sl_int(n_horizontal_pixels);
  pixel_data = new pixel_type[n_rays];
  return pixel_data;
}

void CCPi::instrument::set_pixel_data(pixel_type *p, const sl_int n)
{
  sl_int n_rays = sl_int(n_angles) * sl_int(n_vertical_pixels)
    * sl_int(n_horizontal_pixels);
  if (n != n_rays)
    std::cerr << "Size mismatch setting pixel data\n";
  pixel_data = p;
}

void CCPi::instrument::set_phi(real *p, const int n)
{
  phi = p;
  n_angles = n;
}

bool CCPi::cone_beam::supports_distributed_memory() const
{
  return false;
}

bool CCPi::parallel_beam::supports_distributed_memory() const
{
  return true;
}

bool CCPi::cone_beam::supports_blocks() const
{
  return false;
}

bool CCPi::parallel_beam::supports_blocks() const
{
  return true;
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
    instrument::backward_project(get_h_pixels(), get_v_pixels(), get_phi(),
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
    instrument::backward_project(get_h_pixels(), get_v_pixels(), get_phi(),
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
  for (sl_int curr_ray_z = 0; curr_ray_z < n_rays_z; curr_ray_z++) {
    //real z = det_z[curr_ray_z];
    sl_int k = (sl_int)std::floor((det_z[curr_ray_z]
			       - grid_offset[2]) / voxel_size[2]);
    if (k < 0 or k >= nz_voxels)
      continue;
    sl_int k_offset = k * nx_voxels * ny_voxels;
    sl_int z_offset = curr_ray_z * n_rays_y;

    for (sl_int i = 0; i < matrix_size; i++)
      ray_data[forward_rows[i] + z_offset] += forward_matrix[i]
	* vol_data[k_offset + forward_cols[i]];
    /*
    char desc[6];
    desc[0] = 'G';
    desc[3] = 'C';
    sl_int n_rays = sl_int(n_angles) * sl_int(n_rays_z) * sl_int(n_rays_y);
    sl_int n_vox = sl_int(nx_voxels) * sl_int(ny_voxels) * sl_int(nz_voxels);
    real alpha = 1.0;
    real beta = 1.0;
    sl_int sz = matrix_size;
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
  for (sl_int k = 0; k < nz_voxels; k++) {
    sl_int k_offset = k * nx_voxels * ny_voxels;
    sl_int z_min = -1;
    sl_int z_max = -1;
    // Todo - calculate range rather than looping
    for (sl_int curr_ray_z = 0; curr_ray_z < n_rays_z; curr_ray_z++) {
      sl_int kz = (sl_int)std::floor((det_z[curr_ray_z]
				  - grid_offset[2]) / voxel_size[2]);
      if (kz == k) {
	if (z_min == -1)
	  z_min = curr_ray_z;
	z_max = curr_ray_z;
      } else if (kz > k)
	break;
    }

    sl_int size = nx_voxels * ny_voxels;
    for (sl_int curr_ray_z = z_min; curr_ray_z <= z_max; curr_ray_z++) {
      sl_int z_offset = curr_ray_z * n_rays_y;
      for (sl_int i = 0; i < size; i++) {
	for (sl_int j = backward_rowb[i]; j < backward_rowe[i]; j++) {
	  vol_data[k_offset + i] += backward_matrix[j] *
	    ray_data[backward_cols[j] + z_offset];
	}
      }
      /*
      char desc[6];
      desc[0] = 'G';
      desc[3] = 'C';
      sl_int n_rays = sl_int(n_angles) * sl_int(n_rays_z) * sl_int(n_rays_y);
      real alpha = 1.0;
      real beta = 1.0;
      mkl_dcsrmv("N", &size, &n_rays, &alpha, desc, backward_matrix,
		 backward_cols, backward_rowb, backward_rowe,
		 &ray_data[z_offset], &beta, &vol_data[k_offset]);
      */
    }
  }
}
