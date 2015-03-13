
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

int CCPi::cone_beam::get_z_size(const int n, const real hsize,
				const real vsize) const
{
  int nangles = get_num_angles();
  int nh = get_num_h_pixels();
  int nv = total_num_v_pixels();
  int nz = 0;
  // scan through possibilities to make sure its big enough
  const real_1d &v_range = get_v_pixels();
  const real_1d &h_range = get_h_pixels();
  const real_1d &angles = get_phi();
  real p1z = source_z;
  real p2z = v_range[nv - 1];
  real xy = 0.5 * n * hsize;
  for (int a = 0; a < nangles; a++) {
    real cos_a = std::cos(angles[a]);
    real sin_a = std::sin(angles[a]);
    real p1x = cos_a * source_x - sin_a * source_y;
    real p1y = sin_a * source_x + cos_a * source_y;
    for (int h = 0; h < nh; h++) {
      real p2x = cos_a * detector_x - sin_a * h_range[h];
      real p2y = sin_a * detector_x + cos_a * h_range[h];
      // calc 2D alpha then see what z this would be
      real ax0 = (-xy - p1x) / (p2x - p1x);
      real axn = (xy - p1x) / (p2x - p1x);
      real ay0 = (-xy - p1y) / (p2y - p1y);
      real ayn = (xy - p1y) / (p2y - p1y);
      real ax = std::max(ax0, axn);
      real ay = std::max(ay0, ayn);
      real alpha = std::min(std::min(ax, ay), 1.0);
      real z = alpha * (p2z - p1z) + p1z;
      int m = int(std::ceil(z / vsize) * 2.0);
      if (m > nz)
	nz = m;
    }
  }
  return calc_v_alignment(nz, true);
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
