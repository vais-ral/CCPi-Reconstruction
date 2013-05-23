
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
  return n_angles * n_vertical_pixels * n_horizontal_pixels;
}

pixel_type *const CCPi::instrument::get_pixel_data() const
{
  return pixel_data;
}

pixel_type *CCPi::instrument::create_pixel_data()
{
  long n_rays = n_angles * n_vertical_pixels * n_horizontal_pixels;
  pixel_data = new pixel_type[n_rays];
  return pixel_data;
}

void CCPi::instrument::set_pixel_data(pixel_type *p, const long n)
{
  long n_rays = n_angles * n_vertical_pixels * n_horizontal_pixels;
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
  instrument::backward_project(get_h_pixels(), get_v_pixels(), get_phi(),
			       get_pixel_data(), voxels,
			       get_num_angles(), get_num_h_pixels(),
			       get_num_v_pixels(), origin, width, nx, ny, nz);
  bptime.accumulate();
  bptime.output("backward projection");
}
