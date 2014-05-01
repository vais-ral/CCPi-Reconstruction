
#include <iostream>
#include <omp.h>
#include <cmath>
#ifdef MATLAB_MEX_FILE
#  include "mex_types.hpp"
#else
#  include "base_types.hpp"
#endif // mex
#include "instruments.hpp"
#include "project_line.hpp"
#include "parallel_b.hpp"
#include "parallel_f.hpp"
#include "timer.hpp"

#ifndef USE_TIMER
#  define USE_TIMER false
#endif // USE_TIMER

bool CCPi::parallel_beam::supports_distributed_memory() const
{
  return true;
}

bool CCPi::parallel_beam::supports_blocks() const
{
  return true;
}

void CCPi::parallel_beam::safe_forward_project(pixel_type *pixels,
					       voxel_type *const voxels,
					       const real origin[3],
					       const real width[3],
					       const int nx, const int ny,
					       const int nz) const
{
  timer fptime(USE_TIMER);
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
    /*
    instrument::backward_project(get_h_pixels(), get_v_pixels(), get_phi(),
				 pixels, voxels, get_num_angles(),
				 get_num_h_pixels(), get_num_v_pixels(), origin,
				 width, nx, ny, nz);
    */ /**/
  b2D(get_h_pixels(), get_v_pixels(), get_phi(),
      pixels, voxels, get_num_angles(),
      get_num_h_pixels(), get_num_v_pixels(), origin,
      width, nx, ny, nz);
  /**/
  bptime.accumulate();
  bptime.output("backward projection");
}

void CCPi::parallel_beam::backward_project(voxel_type *const voxels,
					   const real origin[3],
					   const real width[3], const int nx,
					   const int ny, const int nz) const
{
  timer bptime(USE_TIMER);
    /*
    instrument::backward_project(get_h_pixels(), get_v_pixels(), get_phi(),
				 get_pixel_data(), voxels,
				 get_num_angles(), get_num_h_pixels(),
				 get_num_v_pixels(), origin, width, nx, ny, nz);
    */ /**/
  b2D(get_h_pixels(), get_v_pixels(), get_phi(),
      get_pixel_data(), voxels,
      get_num_angles(), get_num_h_pixels(),
      get_num_v_pixels(), origin, width, nx, ny, nz);
  /**/
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
  /**/
  f2D(get_h_pixels(), get_v_pixels(), get_phi(), get_num_angles(),
      get_num_h_pixels(), get_num_v_pixels(), origin, width, nx, ny, nz,
      pixels, voxels);
  /**/ /*
    instrument::forward_project(get_h_pixels(), get_v_pixels(), get_phi(),
				pixels, voxels, get_num_angles(),
				get_num_h_pixels(), get_num_v_pixels(), origin,
				width, nx, ny, nz);
  */
  fptime.accumulate();
  fptime.output(" forward projection");
}
