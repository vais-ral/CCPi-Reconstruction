
#include <iostream>
#include <omp.h>
#include <cmath>
#include <algorithm>
#ifdef MATLAB_MEX_FILE
#  include "mex_types.hpp"
#else
#  include "base_types.hpp"
#endif // mex
#include "instruments.hpp"
#include "project_line.hpp"
#include "cone_b.hpp"
#include "cone_f.hpp"
#include "timer.hpp"

#ifndef USE_TIMER
#  define USE_TIMER false
#endif // USE_TIMER

bool CCPi::cone_beam::supports_distributed_memory() const
{
  return false;
}

bool CCPi::cone_beam::supports_blocks() const
{
  return false;
}

void CCPi::cone_beam::set_params(const real sx, const real sy, const real sz,
				 const real dx, real dy[], real dz[],
				 real ang[], const int ny, const int nz,
				 const int nang)
{
  set_source(sx, sy, sz);
  set_detector(dx);
  real_1d &hp = set_h_pixels(ny);
  for (int i = 0; i < ny; i++)
    hp[i] = dy[i];
  real_1d &vp = set_v_pixels(nz);
  for (int i = 0; i < nz; i++)
    vp[i] = dz[i];
  real_1d &phi = set_phi(nang);
  for (int i = 0; i < nang; i++)
    phi[i] = ang[i];
}

void CCPi::cone_beam::safe_forward_project(pixel_data &pixels,
					   voxel_data &voxels,
					   const real origin[3],
					   const real width[3],
					   const int nx, const int ny,
					   const int nz)
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

void CCPi::cone_beam::forward_project(pixel_data &pixels,
				      voxel_data &voxels,
				      const real origin[3],
				      const real width[3], const int nx,
				      const int ny, const int nz)
{
  timer fptime(USE_TIMER);
  /*
  instrument::forward_project(source_x, source_y, source_z, detector_x,
			      get_h_pixels(), get_v_pixels(), get_phi(),
			      pixels, voxels, get_num_angles(),
			      get_num_h_pixels(), get_num_v_pixels(), origin,
			      width, nx, ny, nz);
  */
  f2D(source_x, source_y, source_z, detector_x, get_h_pixels(), get_v_pixels(),
      get_phi(), pixels, voxels, get_num_angles(), get_num_h_pixels(),
      get_num_v_pixels(), origin, width, nx, ny, nz);
  fptime.accumulate();
  fptime.output(" forward projection");
}

void CCPi::cone_beam::backward_project(pixel_data &pixels,
				       voxel_data &voxels,
				       const real origin[3],
				       const real width[3], const int nx,
				       const int ny, const int nz)
{
  timer bptime(USE_TIMER);
  /*
  instrument::backward_project(source_x, source_y, source_z, detector_x,
			       get_h_pixels(), get_v_pixels(), get_phi(),
			       pixels, voxels, get_num_angles(),
			       get_num_h_pixels(), get_num_v_pixels(), origin,
			       width, nx, ny, nz);
  */
  b2D(source_x, source_y, source_z, detector_x, get_h_pixels(), get_v_pixels(),
      get_phi(), pixels, voxels, get_num_angles(), get_num_h_pixels(),
      get_num_v_pixels(), origin, width, nx, ny, nz);
  bptime.accumulate();
  bptime.output("backward projection");
}

void CCPi::cone_beam::backward_project(voxel_data &voxels,
				       const real origin[3],
				       const real width[3], const int nx,
				       const int ny, const int nz)
{
  timer bptime(USE_TIMER);
  /*
  instrument::backward_project(source_x, source_y, source_z, detector_x,
			       get_h_pixels(), get_v_pixels(), get_phi(),
			       get_pixel_data(), voxels,
			       get_num_angles(), get_num_h_pixels(),
			       get_num_v_pixels(), origin, width, nx, ny, nz);
  */
  b2D(source_x, source_y, source_z, detector_x, get_h_pixels(), get_v_pixels(),
      get_phi(), get_pixel_data(), voxels, get_num_angles(), get_num_h_pixels(),
      get_num_v_pixels(), origin, width, nx, ny, nz);
  bptime.accumulate();
  bptime.output("backward projection");
}
