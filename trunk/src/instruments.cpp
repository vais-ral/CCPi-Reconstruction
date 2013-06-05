
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
		  get_num_h_pixels(), origin, width, nx, ny, nz);
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

    int offset = 0;
    for (long curr_angle = 0; curr_angle < n_angles; curr_angle++) {

      long ray_offset = curr_angle * long(n_rays_y) * long(n_rays_z)
	+ curr_ray_z * long(n_rays_y);

      /* loop over y values on detector */
      for (long curr_ray_y = 0; curr_ray_y < n_rays_y; curr_ray_y++) {

	real data = 0.0;
	for (int i = 0; i < forward_sizes[offset]; i++)
	  data += vol_data[k_offset + forward_y[offset][i] * nx_voxels
			   + forward_x[offset][i]] * forward_data[offset][i];
	ray_data[ray_offset + curr_ray_y] = data;
	offset++;
      }
    }
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
    /*
    // Todo - improve by calculating curr_ray_z range rather than 2 loops
    // we loop on k to avoid possible any possible shared access to voxels
    // Todo - if we stick with 2 loops is there a better loop order for k/z/j/i?
    // just calc k range here then loop over it inside loop over p?
    long k_offset = k * nx_voxels * ny_voxels;
    for (long curr_ray_z = 0; curr_ray_z < n_rays_z; curr_ray_z++) {
      long kz = (long)std::floor((det_z[curr_ray_z]
				  - grid_offset[2]) / voxel_size[2]);
      if (kz != k)
	continue;

      int offset = 0;
      for (long j = 0; j < ny_voxels; j++) {
	long j_offset = k_offset + j * nx_voxels;
	for (long i = 0; i < nx_voxels; i++) {
	  real data = 0.0;
	  for (int p = 0; p < backward_sizes[offset]; p++)
	    data += ray_data[backward_angles[offset][p] * long(n_rays_z)
			     * long(n_rays_y) + curr_ray_z * long(n_rays_y)
			     + backward_h[offset][p]]
	      * backward_data[offset][p];
	  vol_data[j_offset + i] = data;
	  offset++;
	}
      }
    }
    */
    long k_offset = k * nx_voxels * ny_voxels;
    long z_min = -1;
    long z_max = -1;
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

    int offset = 0;
    for (long j = 0; j < ny_voxels; j++) {
      long j_offset = k_offset + j * nx_voxels;
      for (long i = 0; i < nx_voxels; i++) {
	real data = 0.0;
	for (int p = 0; p < backward_sizes[offset]; p++) {
	  long step = backward_angles[offset][p] * long(n_rays_z)
	    * long(n_rays_y) + z_min * long(n_rays_y) + backward_h[offset][p];
	  for (long curr_ray_z = z_min; curr_ray_z <= z_max; curr_ray_z++) {
	    data += ray_data[step] * backward_data[offset][p];
	    step += long(n_rays_y);
	  }
	}
	vol_data[j_offset + i] = data;
	offset++;
      }
    }
  }
}
