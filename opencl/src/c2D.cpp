
#include <float.h>
#ifdef MATLAB_MEX_FILE
#  include "mex_types.hpp"
#else
#  include "base_types.hpp"
#endif // mex
#include "instruments.hpp"
#include "timer.hpp"
#include "ui_calls.hpp"
#include "accel.hpp"
#if defined(TEST2D) || defined(CBZCHECK) 
#  include <iostream>
#endif // TEST2D
#include <omp.h>

static const recon_type epsilon = FLT_EPSILON;

void CCPi::cone_beam::calc_xy_z(pixel_type *const pixels,
				const voxel_ptr_1d &voxels,
				const recon_1d &alpha_xy, const int n,
				const recon_type pzbz, const int nv,
				const int nz, const int midp,
				const recon_1d &delta_z,
				const recon_1d &inv_delz, const recon_1d &vox_z,
				int_1d &kv)
{
  // alpha always increases
  // delta_z[i] = v_pixels[i] - source_z
  // inv_delz[i] = 1.0 / delta_z[i] = 1.0 / (v_pixels[i] - source_z)
  // vox_z[i] = vox_origin[2] + real(i) * vox_size[2] - source_z
  // so alpha_z is intercept of line at voxel[k] position
  // (vox_origin[2] + k * vox_size[2]) = p1_z + alpha_z * (v_pix[v] - p1_z)
  // alpha_z = ((vox_origin[2] + k * vox_size[2]) - p1_z) / (v_pix[v] - p1_z)
  // alpha_z = vox_z[k] / delta_z[v] = vox_z[k] * inv_delz[v]
  // and for k calc
  // k = int((p1_z + alpha_xy[m - 1] * (v_pix[v] - p1_z) - b_z) / d_z)
  // k = int((p1_z + alpha_xy[m - 1] * (p2_z[v] - p1_z) - b_z) / d_z)
  // k = int((p1_z + alpha_xy[m - 1] * delta_z[v] - b_z) * inv_dz)
  // k = int((p1_z - b_z) * inv_dz + (alpha_xy[m - 1] * delta_z[v]) * inv_dz)
  // k = int((pzbz + (alpha_xy[m - 1] * inv_dz) * delta_z[v])
  // k = int((pzbz + alpha_inv * delta_z[v])
  recon_type *dz_ptr = assume_aligned(&(delta_z[0]), recon_type);
  recon_type *iz_ptr = assume_aligned(&(inv_delz[0]), recon_type);
  recon_type *vz_ptr = assume_aligned(&(vox_z[0]), recon_type);
  recon_type *axy_ptr = assume_aligned(&(alpha_xy[0]), recon_type);
  int *k_ptr = assume_aligned(&(kv[0]), int);
  pixel_type *const pix = assume_aligned(pixels, pixel_type);
  recon_type alpha_m0 = axy_ptr[0];
  for (int m = 1; m < n; m++) {
    const voxel_type *const vox = assume_aligned(voxels[m], voxel_type);
    for (int v = 0; v < nv; v++)
      k_ptr[v] = int(pzbz + alpha_m0 * dz_ptr[v]);
    const recon_type alpha_m1 = axy_ptr[m];
    for (int v = 0; v < midp; v++) {
      int k = k_ptr[v];
      recon_type alpha_z = vz_ptr[k] * iz_ptr[v];
      recon_type min_z = std::min(alpha_z, alpha_m1);
      pix[v] += (vox[k] * (min_z - alpha_m0)
		 + vox[k - 1] * (alpha_m1 - min_z));
    }
    for (int v = midp; v < nv; v++) {
      int k = k_ptr[v];
      recon_type alpha_z = vz_ptr[k + 1] * iz_ptr[v];
      recon_type min_z = std::min(alpha_z, alpha_m1);
      pix[v] += (vox[k] * (min_z - alpha_m0)
		 + vox[k + 1] * (alpha_m1 - min_z));
    }
    alpha_m0 = alpha_m1;
  }
}

void CCPi::cone_beam::fproject_xy(const real p1_x, const real p1_y,
				  const real p2_x, const real p2_y,
				  pixel_data &pixels, voxel_data &voxels,
				  const real b_x, const real b_y,
				  const real d_x, const real d_y, const int nx,
				  const int ny, const int nz, const int a,
				  const int h, const int nv,
				  const sl_int ij_base, const sl_int nyz,
				  recon_1d &alpha_xy, voxel_ptr_1d &ij_arr,
				  int_1d &ij_index, const sl_int block_yz,
				  const int xy_base, int &count)
{
  // Find intercepts with voxels (x1,y1) (x2, y2)
  // line goes from (p1_x, p1_y) to (p2_x, p2_y) delta = p2 - p1
  // so parametric form x = p1_x + alpha * delta_x, y = p1_y + alpha * delta_y
  // alpha from 0 to 1 contains the voxels.
  // Rescaled so voxel origin at 0, voxels size 1 (px, py) -> (qx, qy)
  // Alpha should be unchanged by scaling since its still the same fraction
  // of the line (I think)
  const real q1_x = (p1_x - b_x) / d_x;
  const real q1_y = (p1_y - b_y) / d_y;
  const real q2_x = (p2_x - b_x) / d_x;
  const real q2_y = (p2_y - b_y) / d_y;
  const real delta_x = q2_x - q1_x;
  const real delta_y = q2_y - q1_y;
  const real inv_dx = 1.0 / delta_x;
  const real inv_dy = 1.0 / delta_y;

  count = 0;
  if (std::abs(delta_x) < epsilon) {
    if (std::abs(delta_y) < epsilon) {
      // Its not a line - shouldn't happen
      //return;
    } else {
      // line parallel to y
      int i = int(std::floor(q1_x)); // == q2_x
      if (i >= 0 and i < nx) {
	if (delta_y < 0.0) {
	  alpha_xy[count] = (real(ny) - q1_y) * inv_dy;
	  voxel_ptr offset = voxels.data() + ij_base + sl_int(i) * nyz
	    + sl_int(ny - 1) * sl_int(nz);
	  int index = xy_base + i * block_yz + (ny - 1) * nz;
	  ij_arr[count] = offset;
	  ij_index[count] = index;
	  count++;
	  for (int j = ny - 1; j >= 0; j--) {
	    alpha_xy[count] = (real(j) - q1_y) * inv_dy;
	    ij_arr[count] = offset;
	    ij_index[count] = index;
	    offset -= nz;
	    index -= nz;
	    count++;
	  }
	} else {
	  alpha_xy[count] = (real(0) - q1_y) * inv_dy;
	  voxel_ptr offset = voxels.data() + ij_base + sl_int(i) * nyz;
	  int index = xy_base + i * block_yz;
	  ij_arr[count] = offset;
	  ij_index[count] = index;
	  count++;
	  for (int j = 0; j < ny; j++) {
	    alpha_xy[count] = (real(j + 1) - q1_y) * inv_dy;
	    ij_arr[count] = offset;
	    ij_index[count] = index;
	    offset += nz;
	    index += nz;
	    count++;
	  }
	}
      }
    }
  } else if (std::abs(delta_y) < epsilon) {
    // line parallel to x
    int j = int(std::floor(q1_y)); // == q2_y
    if (j >= 0 and j < ny) {
      if (delta_x < 0.0) {
	alpha_xy[count] = (real(nx) - q1_x) * inv_dx;
	voxel_ptr offset = voxels.data() + ij_base + sl_int(nx - 1) * nyz
	  + sl_int(j) * sl_int(nz);
	int index = xy_base + (nx - 1) * block_yz;
	ij_arr[count] = offset;
	ij_index[count] = index;
	count++;
	for (int i = nx - 1; i >= 0; i--) {
	  alpha_xy[count] = (real(i) - q1_x) * inv_dx;
	  ij_arr[count] = offset;
	  ij_index[count] = index;
	  offset -= nyz;
	  index -= block_yz;
	  count++;
	}
      } else {
	alpha_xy[count] = (real(0) - q1_x) * inv_dx;
	voxel_ptr offset = voxels.data() + ij_base + sl_int(j) * sl_int(nz);
	int index = xy_base + j * nz;
	ij_arr[count] = offset;
	ij_index[count] = index;
	count++;
	for (int i = 0; i < nx; i++) {
	  alpha_xy[count] = (real(i + 1) - q1_x) * inv_dx;
	  ij_arr[count] = offset;
	  ij_index[count] = index;
	  offset += nyz;
	  index += block_yz;
	  count++;
	}
      }
    }
  } else {
    // general line drawing - find intercepts
    // rescaled so
    const real x_0 = 0.0;
    const real y_0 = 0.0;
    const real x_n = real(nx);
    const real y_n = real(ny);
    const real alpha_x_0 = (x_0 - q1_x) * inv_dx;
    const real alpha_y_0 = (y_0 - q1_y) * inv_dy;
    const real alpha_x_n = (x_n - q1_x) * inv_dx;
    const real alpha_y_n = (y_n - q1_y) * inv_dy;
    const real alpha_x_min = std::min(alpha_x_0, alpha_x_n);
    const real alpha_x_max = std::max(alpha_x_0, alpha_x_n);
    const real alpha_y_min = std::min(alpha_y_0, alpha_y_n);
    const real alpha_y_max = std::max(alpha_y_0, alpha_y_n);
    const real alpha_min = std::max(std::max(alpha_x_min, alpha_y_min),
				    real(0.0));
    const real alpha_max = std::min(std::min(alpha_x_max, alpha_y_max),
				    real(1.0));
    
    if (alpha_min < alpha_max - epsilon) {
      real_1d alpha_x(nx + 1);
      for (int i = 0; i <= nx; i++)
	alpha_x[i] = ((real(i) - q1_x) * inv_dx);
      real_1d alpha_y(ny + 1);
      for (int i = 0; i <= ny; i++)
	alpha_y[i] = ((real(i) - q1_y) * inv_dy);
      int x = 0;
      int y = 0;
      if (delta_x > 0.0) {
	if (delta_y > 0.0) {
	  if (alpha_min == alpha_x_0) {
	    x = 0;
	    y = int(std::floor(q1_y + alpha_min * delta_y));
	  } else if (alpha_min == alpha_y_0) {
	    x = int(std::floor(q1_x + alpha_min * delta_x));
	    y = 0;
	  } else
	    report_error("something wrong in x+ y+");
	  alpha_xy[count] = alpha_min;
	  voxel_ptr offset = voxels.data() + ij_base + sl_int(x) * nyz
	    + sl_int(y) * sl_int(nz);
	  int index = xy_base + x * block_yz + y * nz;
	  ij_arr[count] = offset;
	  ij_index[count] = index;
	  count++;
	  // could do x_next/y_next here and only calc the one that changes
	  // inside the if statements below, which would reduce the flops
	  while (x < nx and y < ny) {
	    // could calc a general alpha_step and maybe have rounding issues
	    ij_arr[count] = offset;
	    ij_index[count] = index;
	    if (alpha_x[x + 1] < alpha_y[y + 1] - epsilon) {
	      alpha_xy[count] = alpha_x[x + 1];
	      offset += nyz;
	      index += block_yz;
	      x++;
	    } else if (alpha_y[y + 1] < alpha_x[x + 1] - epsilon) {
	      alpha_xy[count] = alpha_y[y + 1];
	      offset += nz;
	      index += nz;
	      y++;
	    } else {
	      alpha_xy[count] = std::max(alpha_x[x + 1], alpha_y[y + 1]);
	      offset += nyz + nz;
	      index += block_yz + nz;
	      x++;
	      y++;
	    }
	    count++;
	  }
	} else {
	  if (alpha_min == alpha_y_n) {
	    x = int(std::floor(q1_x + alpha_min * delta_x));
	    y = ny - 1;
	  } else if (alpha_min == alpha_x_0) {
	    x = 0;
	    y = int(std::floor(q1_y + alpha_min * delta_y));
	  } else
	    report_error("something wrong in x+ y-");
	  alpha_xy[count] = alpha_min;
	  voxel_ptr offset = voxels.data() + ij_base + sl_int(x) * nyz
	    + sl_int(y) * sl_int(nz);
	  int index = xy_base + x * block_yz + y * nz;
	  ij_arr[count] = offset;
	  ij_index[count] = index;
	  count++;
	  while (x < nx and y >= 0) {
	    ij_arr[count] = offset;
	    ij_index[count] = index;
	    if (alpha_x[x + 1] < alpha_y[y] - epsilon) {
	      alpha_xy[count] = alpha_x[x + 1];
	      offset += nyz;
	      index += block_yz;
	      x++;
	    } else if (alpha_y[y] < alpha_x[x + 1] - epsilon) {
	      alpha_xy[count] = alpha_y[y];
	      offset -= nz;
	      index -= nz;
	      y--;
	    } else {
	      alpha_xy[count] = std::max(alpha_x[x + 1], alpha_y[y]);
	      offset += nyz - nz;
	      index += block_yz - nz;
	      x++;
	      y--;
	    }
	    count++;
	  }
	}
      } else {
	if (delta_y > 0.0) {
	  if (alpha_min == alpha_x_n) {
	    x = nx - 1;
	    y = int(std::floor(q1_y + alpha_min * delta_y));
	  } else if (alpha_min == alpha_y_0) {
	    x = int(std::floor(q1_x + alpha_min * delta_x));
	    y = 0;
	  } else
	    report_error("something wrong in x- y+");
	  alpha_xy[count] = alpha_min;
	  voxel_ptr offset = voxels.data() + ij_base + sl_int(x) * nyz
	    + sl_int(y) * sl_int(nz);
	  int index = xy_base + x * block_yz + y * nz;
	  ij_arr[count] = offset;
	  ij_index[count] = index;
	  count++;
	  while (x >= 0 and y < ny) {
	    ij_arr[count] = offset;
	    ij_index[count] = index;
	    if (alpha_x[x] < alpha_y[y + 1] - epsilon) {
	      alpha_xy[count] = alpha_x[x];
	      offset -= nyz;
	      index -= block_yz;
	      x--;
	    } else if (alpha_y[y + 1] < alpha_x[x] - epsilon) {
	      alpha_xy[count] = alpha_y[y + 1];
	      offset += nz;
	      index += nz;
	      y++;
	    } else {
	      alpha_xy[count] = std::max(alpha_y[y + 1], alpha_x[x]);
	      offset += -nyz + nz;
	      index += -block_yz + nz;
	      x--;
	      y++;
	    }
	    count++;
	  }
	} else {
	  if (alpha_min == alpha_x_n) {
	    x = nx - 1;
	    if (alpha_min == alpha_y_n)
	      y = ny - 1;
	    else
	      y = int(std::floor(q1_y + alpha_min * delta_y));
	  } else if (alpha_min == alpha_y_n) {
	    x = int(std::floor(q1_x + alpha_min * delta_x));
	    y = ny - 1;
	  } else
	    report_error("something wrong in x- y-");
	  alpha_xy[count] = alpha_min;
	  voxel_ptr offset = voxels.data() + ij_base + sl_int(x) * nyz
	    + sl_int(y) * sl_int(nz);
	  int index = xy_base + x * block_yz + y * nz;
	  ij_arr[count] = offset;
	  ij_index[count] = index;
	  count++;
	  while (x >= 0 and y >= 0) {
	    ij_arr[count] = offset;
	    ij_index[count] = index;
	    if (alpha_x[x] < alpha_y[y] - epsilon) {
	      alpha_xy[count] = alpha_x[x];
	      offset -= nyz;
	      index -= block_yz;
	      x--;
	    } else if (alpha_y[y] < alpha_x[x] - epsilon) {
	      alpha_xy[count] = alpha_y[y];
	      offset -= nz;
	      index -= nz;
	      y--;
	    } else {
	      alpha_xy[count] = std::max(alpha_x[x], alpha_y[y]);
	      offset -= (nyz + nz);
	      index -= (block_yz + nz);
	      x--;
	      y--;
	    }
	    count++;
	  }
	}
      }
      // Todo - check final value is alpha_max?
    }
  }
}

void CCPi::cone_beam::f2D_cpu(const real source_x, const real source_y,
			      const real source_z, const real detector_x,
			      const real_1d &h_pixels, const real_1d &v_pixels,
			      const real_1d &angles, pixel_data &pixels,
			      voxel_data &voxels, const int n_angles,
			      const int n_h, const int n_v,
			      const real grid_offset[3],
			      const real voxel_size[3], const int nx_voxels,
			      const int ny_voxels, const int nz_voxels)
{
  // Todo - check that source_z to det_z max z angle < 45 for safe usage.
  int mid = -1;
  for (int i = 0; i < n_v; i++) {
    if (v_pixels[i] > source_z) {
      mid = i;
      break;
    }
  }
  if (v_pixels[mid - 1] > source_z - epsilon and
      v_pixels[mid - 1] < source_z + epsilon)
    report_error("Oops");
  if (v_pixels[mid] > source_z - epsilon and
      v_pixels[mid] < source_z + epsilon)
    report_error("Oops");

#pragma omp parallel shared(h_pixels, v_pixels, pixels, voxels, angles, voxel_size, grid_offset) firstprivate(n_angles, n_h, n_v, nx_voxels, ny_voxels, nz_voxels, mid)
  {
    int thread_id = omp_get_thread_num();
    int nthreads = omp_get_num_threads();
    // path length from source to detector is independent of rotation
    recon_2d d_conv(boost::extents[n_h][n_v]);
    real x2 = (detector_x - source_x) * (detector_x - source_x);
    for (int i = 0; i < n_h; i++) {
      real xy2 = x2 + (h_pixels[i] - source_y) * (h_pixels[i] - source_y);
      for (int j = 0; j < n_v; j++) {
	d_conv[i][j] = std::sqrt(xy2 + (v_pixels[j] - source_z)
				 * (v_pixels[j] - source_z));
      }
    }

    //const recon_type inv_dz = recon_type(real(1.0) / voxel_size[2]);
    const recon_type pzbz = recon_type((source_z
					- grid_offset[2]) / voxel_size[2]);
    recon_1d delta_z(n_v);
    for (int i = 0; i < n_v; i++)
      delta_z[i] = v_pixels[i] - source_z;
    recon_1d inv_delz(n_v);
    for (int i = 0; i < n_v; i++)
      inv_delz[i] = 1.0 / delta_z[i];
    for (int i = 0; i < n_v; i++)
      delta_z[i] /= voxel_size[2];
    recon_1d vox_z(nz_voxels + 1);
    for (int i = 0; i <= nz_voxels; i++)
      vox_z[i] = grid_offset[2] + real(i) * voxel_size[2] - source_z;

    const real pixel_step = (h_pixels[1] - h_pixels[0]);
    const real ipix_step = 1.0 / pixel_step;
    const real hpix0 = (source_y - h_pixels[0]) / pixel_step;
    const real l = detector_x - source_x;
    sl_int nyz = sl_int(ny_voxels) * sl_int(nz_voxels);
    int max_n = std::max(nx_voxels, ny_voxels);
    recon_1d l_xy(2 * max_n + 2);
    voxel_ptr_1d ij_arr(2 * max_n + 2);
    int_1d ij_offsets(2 * max_n + 2);
    //sl_int nyz = sl_int(ny) * sl_int(nz);
    int nwork = 0;

    const int a_block = n_angles;
    const int x_block = nx_voxels;
    const int y_block = ny_voxels;
    long counter = 0;
    for (int block_a = 0; block_a < n_angles; block_a += a_block) {
      int a_step = a_block;
      if (block_a + a_step > n_angles)
	a_step = n_angles - block_a;
      // we don't block h since we have a min/max from the x/y blocks
      for (int block_x = 0; block_x < nx_voxels; block_x += x_block) {
	int x_step = x_block;
	if (block_x + x_step > nx_voxels)
	  x_step = nx_voxels - block_x;
	sl_int i_base = sl_int(block_x) * nyz;
	real vx = grid_offset[0] + real(block_x) * voxel_size[0];
	real wx = grid_offset[0] + real(block_x + x_step) * voxel_size[0];
	for (int block_y = 0; block_y < ny_voxels; block_y += y_block) {
	  int y_step = y_block;
	  if (block_y + y_step > ny_voxels)
	    y_step = ny_voxels - block_y;
	  sl_int ij_base = i_base + sl_int(block_y) * sl_int(nz_voxels);
	  real vy = grid_offset[1] + real(block_y) * voxel_size[1];
	  real wy = grid_offset[1] + real(block_y + y_step) * voxel_size[1];
	  for (int ax = 0; ax < a_step; ax++) {
	    int max_n = std::max(nx_voxels, ny_voxels);
	    int_1d kv(n_v);
	    int a = block_a + ax;
	    real cphi = std::cos(angles[a]);
	    real sphi = std::sin(angles[a]);
	    real p1_x = cphi * source_x - sphi * source_y;
	    real p1_y = sphi * source_x + cphi * source_y;

	    // from bproj
	    //int hmin = 0;
	    //int hmax = n_h - 1;
	    /**/
	    real y00;
	    real y01;
	    real y10;
	    real y11;
	    const real delta_x_0 = vx - p1_x;
	    const real delta_y_0 = vy - p1_y;
	    const real delta_x_n = wx - p1_x;
	    const real delta_y_n = wy - p1_y;
	    const real cx_0 = cphi * vx - source_x; // Todo - could pass this in
	    const real cx_n = cphi * wx - source_x; // cx_0 + cphi * d_x
	    const real sy_0 = sphi * vy; // Todo could pass in as array
	    const real sy_n = sphi * wy;
	    if (std::abs(sphi) > std::abs(cphi)) {
	      const real ilphi = l / sphi;
	      real u = sy_0 + cx_0;
	      y00 = (cphi - delta_x_0 / u) * ilphi;
	      u = sy_n + cx_0;
	      y01 = (cphi - delta_x_0 / u) * ilphi;
	      u = sy_0 + cx_n;
	      y10 = (cphi - delta_x_n / u) * ilphi;
	      u = sy_n + cx_n;
	      y11 = (cphi - delta_x_n / u) * ilphi;
	    } else {
	      real ilphi = l / cphi;
	      real u = cx_0 + sy_0;
	      y00 = (-sphi + delta_y_0 / u) * ilphi;
	      u = cx_0 + sy_n;
	      y01 = (-sphi + delta_y_n / u) * ilphi;
	      u = cx_n + sy_0;
	      y10 = (-sphi + delta_y_0 / u) * ilphi;
	      u = cx_n + sy_n;
	      y11 = (-sphi + delta_y_n / u) * ilphi;
	    }
	    int h00 = int(std::floor(y00 * ipix_step + hpix0));
	    int h01 = int(std::floor(y01 * ipix_step + hpix0));
	    int hmin = std::min(h00, h01);
	    int hmax = std::max(h00, h01);
	    int h10 = int(std::floor(y10 * ipix_step + hpix0));
	    hmin = std::min(hmin, h10);
	    hmax = std::max(hmax, h10);
	    int h11 = int(std::floor(y11 * ipix_step + hpix0));
	    hmin = std::max(std::min(hmin, h11), 0);
	    hmax = std::min(std::max(hmax, h11), n_h - 1);
	    /**/
	    // end bproj
	  
	    for (int h = hmin; h <= hmax; h++) {
	      if (counter % nthreads == thread_id) {
		// rotate source and detector positions by current angle
		real p2_x = cphi * detector_x - sphi * h_pixels[h];
		real p2_y = sphi * detector_x + cphi * h_pixels[h];
		fproject_xy(p1_x, p1_y, p2_x, p2_y, pixels, voxels, vx, vy,
			    voxel_size[0], voxel_size[1], x_step, y_step,
			    nz_voxels, a, h, n_v, ij_base, nyz, l_xy, ij_arr,
			    ij_offsets, nyz, ij_base, nwork);
		if (nwork > 0) {
		  if (nwork > 2 * max_n + 2)
		    report_error("forward project overflow", nwork);
		  else {
#ifdef CBZCHECK
		    const int nzm1 = nz_voxels - 1;
		    int min_xy_all = nwork;
		    for (int m = 1; m < nwork; m++) {
		      int k = int(std::floor(pzbz + l_xy[m - 1] * delta_z[0]));
		      if (k <= 0) {
			if (k == 0)
			  min_xy_all = m;
			else {
			  min_xy_all = 0;
			}
			break;
		      }
		    }
		    if (min_xy_all < nwork)
		      std::cerr << "Upper z extent wrong for cone-beam "
				<< min_xy_all << ' ' << nwork << '\n';
		    int max_xy_all = nwork;
		    for (int m = 1; m < nwork; m++) {
		      int k = int(std::floor(pzbz + l_xy[m - 1]
					     * delta_z[n_v - 1]));
		      if (k >= nzm1) {
			if (k == nzm1)
			  max_xy_all = m;
			else {
			  max_xy_all = 0;
			}
			break;
		      }
		    }
		    if (max_xy_all < nwork)
		      std::cerr << "Upper z extent wrong for cone-beam "
				<< max_xy_all << ' ' << nwork << '\n';
#endif // CBZCHECK
		    calc_xy_z(&(pixels[a][h][0]), ij_arr, l_xy, nwork, pzbz,
			      n_v, nz_voxels, mid, delta_z, inv_delz, vox_z,
			      kv);
		  }
		}
	      }
	      counter++;
	    }
	  }
	}
      }
      // Todo - check that this is ok, 0.58 -> 0.63 in vis
      int a_idx = a_step / nthreads;
      if (a_step % nthreads != 0)
	a_idx++;
      int a_off = a_idx * thread_id;
      if (a_off + a_idx > a_step)
	a_idx = a_step - a_off;
      for (int ax = 0; ax < a_idx; ax++) {
	int a = block_a + a_off + ax;
	const recon_type *const dc = assume_aligned(&(d_conv[0][0]),
						    recon_type);
	pixel_type *const pix = assume_aligned(&(pixels[a][0][0]), pixel_type);
	for (int hv = 0; hv < n_h * n_v; hv++)
	  pix[hv] *= dc[hv];
      }
    }
  }
}

void CCPi::cone_beam::f2D(const real source_x, const real source_y,
			  const real source_z, const real detector_x,
			  const real_1d &h_pixels, const real_1d &v_pixels,
			  const real_1d &angles, pixel_data &pixels,
			  voxel_data &voxels, const int n_angles,
			  const int n_h, const int n_v,
			  const real grid_offset[3], const real voxel_size[3],
			  const int nx_voxels, const int ny_voxels,
			  const int nz_voxels)
{
  if (machine::has_accelerator())
    f2D_accel(source_x, source_y, source_z, detector_x, h_pixels, v_pixels,
	      angles, pixels, voxels, n_angles, n_h, n_v, grid_offset,
	      voxel_size, nx_voxels, ny_voxels, nz_voxels);
  else
    f2D_cpu(source_x, source_y, source_z, detector_x, h_pixels, v_pixels,
	    angles, pixels, voxels, n_angles, n_h, n_v, grid_offset,
	    voxel_size, nx_voxels, ny_voxels, nz_voxels);          
}

/*
  vox_z[k] = b_z + k * d_z
  v_pix[v] = v_pix[0] + v * v_step -> v_step = v_pix[1] - v_pix[0]
  so for a parametric line with z = source_z + alpha * (v_pix[v] - source_z)
  we can solve for a particular vox_z[k] by
    vox_z[k] = source_z + alpha * (v_pix[v] - source_z)
  which ends up as
    v = (source_z - v_pix[0]) + (vox_z[k] - source_z)
            -------------        -------------------
              v_step               alpha * v_step
*/

// Todo - loop over k and find v? alpha_xy may be the wrong info for this.
void CCPi::cone_beam::calc_ah_z(const pixel_ptr_1d &pixels,
				voxel_type *const voxels,
				const recon_1d &alpha_xy_0,
				const recon_1d &alpha_xy_1, const int n,
				const recon_type pzbz, const int nv,
				const int nz, const int midp,
				const recon_1d &delta_z,
				const recon_1d &inv_delz, const recon_1d &vox_z,
				int_1d &kv)
{
  // pzdv = (p1z - vpix[0]) / vpix_step
  // z_1 = (bz + 1 * dz - p1z) / vpix_step ,
  // z_nm = (bz + (nz - 1) * dz - p1z) / vpix_step
  // We have rounding issues with pzdv, so add a small increment
  // since epsilon is 1 + epsilon, we need to scale so its meaningful
  recon_type *dz_ptr = assume_aligned(&(delta_z[0]), recon_type);
  recon_type *iz_ptr = assume_aligned(&(inv_delz[0]), recon_type);
  recon_type *vz_ptr = assume_aligned(&(vox_z[0]), recon_type);
  recon_type *axy0_ptr = assume_aligned(&(alpha_xy_0[0]), recon_type);
  recon_type *axy1_ptr = assume_aligned(&(alpha_xy_1[0]), recon_type);
  int *k_ptr = assume_aligned(&kv[0], int);
  voxel_type *const vox = assume_aligned(voxels, voxel_type);
  for (int m = 0; m < n; m++) {
    const pixel_type *const pix = assume_aligned(pixels[m], pixel_type);
    const recon_type alpha_m0 = axy0_ptr[m];
    const recon_type alpha_m1 = axy1_ptr[m];
    for (int v = 0; v < nv; v++)
      k_ptr[v] = int(pzbz + alpha_m0 * dz_ptr[v]);
    for (int v = 0; v < midp; v++) {
      int k = k_ptr[v];
      recon_type alpha_z = vz_ptr[k] * iz_ptr[v];
      recon_type min_z = std::min(alpha_z, alpha_m1);
      vox[k] += pix[v] * (min_z - alpha_m0);
      vox[k - 1] += pix[v] * (alpha_m1 - min_z);
    }
    for (int v = midp; v < nv; v++) {
      int k = k_ptr[v];
      recon_type alpha_z = vz_ptr[k + 1] * iz_ptr[v];
      recon_type min_z = std::min(alpha_z, alpha_m1);
      vox[k] += pix[v] * (min_z - alpha_m0);
      vox[k + 1] += pix[v] * (alpha_m1 - min_z);
    }
  }
}

/*
  At no rotation the line from the source at the centre of the cone
  is (sx, sy) to (dx, sy). In rotated form this is
      (x) = (sx * cphi - sy * sphi) + t ((dx - sx) cphi)
      (y) = (sx * sphi + sy * cphi) +   ((dx - sx) sphi)
  rescaling t
      (x) = (sx * cphi - sy * sphi) + t (cphi)
      (y) = (sx * sphi + sy * cphi) +   (sphi)
  so the line representing the detector pixels which is perpendicular to this is
      (x) = (dx * cphi - sy * sphi) + t ( sphi)
      (y) = (dx * sphi + sy * cphi) +   (-cphi)
  Let the source rotated though theta be (px, py) , then a line from the
  detector through a voxel point (vx, vy) (such as its centre) is
      (x) = (px) + u (vx - px)
      (y) = (py) +   (vy - py)
  If we use (q_x, q_y) for the initial point in the detector line, then
  this intercepts the detector line at
      (px) + u (vx - px) = (qx) + t ( sphi)
      (py) +   (vy - py) = (qy) +   (-cphi)
  where u is not equal to 0 since that would be the source position.
  If cphi != 0 then
      t = (qy - u (vy-py) - py) / cphi
    and
      px + u (vx-px) = qx + (qy - u (vy-py) - py) sphi / cphi
    so
      u [(vx-px) + (vy-py) sphi/cphi] = qx - px + (qy - py) sphi/cphi
  Alternatively if sphi != 0
      t = (px + u (vx-px) - qx) / sphi
      py + u (vy-px) = qy - (px + u (vx-px) - qx) cphi / sphi
   => u [(vy-py) + (vx-px) cphi/sphi] = qy - py + (qx - px) cphi/sphi
  So we solve for t, generate u and apply u.
  Then the detector line transformed back to 0 angle has
  x = dx and y = sy - t
*/

void CCPi::cone_beam::bproject_ah(const real source_x, const real source_y,
				  pixel_data &pixels, voxel_data &voxels,
				  const real x_0, const real y_0,
				  const real x_n, const real y_n, const int nz,
				  const int i, const int j, const int n_angles,
				  const int n_h, const int n_v,
				  const real_1d &h_pixels, const int midp,
				  const real_1d &cangle, const real_1d &sangle,
				  const recon_1d &delta_z,
				  const recon_1d &inv_delz,
				  const recon_1d &vox_z, const recon_type pzbz,
				  const recon_type pzdv, const recon_type z_1,
				  const recon_type z_nm, const real_1d &p1x,
				  const real_1d &p1y, const real_1d &cdetx,
				  const real_1d &sdetx, const real_1d &ilcphi,
				  const real_1d &ilsphi, const int a_off,
				  pixel_ptr_1d &ah_arr, recon_1d &alpha_xy_0,
				  recon_1d &alpha_xy_1, int_1d &ah_index,
				  const int a_base, int &count)
{
  // Rather than using the centre just calculate for all 4 corners,
  // generate h values and loop from smallest to largest.
  count = 0;
  // corners (x0,y0), (x0,yn), (xn,y0), (xn,yn)
  const real pixel_step = h_pixels[1] - h_pixels[0];
  const real ipix_step = 1.0 / pixel_step;
  const real hpix0 = (source_y - h_pixels[0]) / pixel_step;
  sl_int nah = sl_int(n_h) * sl_int(n_v);
  pixel_ptr ah_offset = pixels.data() + sl_int(a_off) * nah;
  int ah_pos = a_base * nah;
  for (int ax = 0; ax < n_angles; ax++) {
    int a = a_off + ax;
    const real cphi = cangle[a];
    const real sphi = sangle[a];
    const real p1_x = p1x[a];
    const real p1_y = p1y[a];
    //const real qdx = cdetx[a];
    //const real qdy = sdetx[a];
    //const real qx = qdx - spy[a];
    //const real qy = qdy + cpy[a];
    real y00;
    real y01;
    real y10;
    real y11;
    //const real qpx = qx - p1_x;
    //const real qpy = qy - p1_y;
    const real delta_x_0 = x_0 - p1_x;
    const real delta_y_0 = y_0 - p1_y;
    const real delta_x_n = x_n - p1_x;
    const real delta_y_n = y_n - p1_y;
    const real cx_0 = cphi * x_0 - source_x; // Todo - could pass this in
    const real cx_n = cphi * x_n - source_x; // cx_0 + cphi * d_x
    const real sy_0 = sphi * y_0; // Todo could pass in as array
    const real sy_n = sphi * y_n;
    if (std::abs(sphi) > std::abs(cphi)) {
      const real ilphi = ilsphi[a];
      //real up = l / (sphi * y_0 + cphi * x_0 - source_x);
      real u = sy_0 + cx_0;
      //real tp = (- cphi * l + delta_x_0 * u) /sphi;
      //real t = (- cphi + delta_x_0 / u) * ilphi;
      //y00 = source_y - t;
      y00 = (cphi - delta_x_0 / u) * ilphi;
      //up = l / (sphi * y_n + cphi * x_0 - source_x);
      u = sy_n + cx_0;
      //tp = (- cphi * l + delta_x_0 * u) / sphi;
      //t = (- cphi + delta_x_0 / u) * ilphi;
      //y01 = source_y - t;
      y01 = (cphi - delta_x_0 / u) * ilphi;
      //up = l / (sphi * y_0 + cphi * x_n - source_x);
      u = sy_0 + cx_n;
      //tp = (- cphi * l + delta_x_n * u) / sphi;
      //t = (- cphi + delta_x_n / u) * ilphi;
      //y10 = source_y - t;
      y10 = (cphi - delta_x_n / u) * ilphi;
      //up = l / (sphi * y_n + cphi * x_n - source_x);
      u = sy_n + cx_n;
      //tp = (- cphi * l + delta_x_n * u) / sphi;
      //t = (- cphi + delta_x_n / u) * ilphi;
      //y11 = source_y - t;
      y11 = (cphi - delta_x_n / u) * ilphi;
    } else {
      real ilphi = ilcphi[a];
      //real up = cphi * x_0 + sphi * y_0 - source_x;
      real u = cx_0 + sy_0;
      //real tp = (sphi - delta_y_0 / up) * ilphi;
      //real t = (sphi - delta_y_0 / u) * ilphi;
      //y00 = source_y - t;
      y00 = (-sphi + delta_y_0 / u) * ilphi;
      u = cx_0 + sy_n;
      //t = (sphi - delta_y_n / u) * ilphi;
      //y01 = source_y - t;
      y01 = (-sphi + delta_y_n / u) * ilphi;
      u = cx_n + sy_0;
      //t = (sphi - delta_y_0 / u) * ilphi;
      //y10 = source_y - t;
      y10 = (-sphi + delta_y_0 / u) * ilphi;
      u = cx_n + sy_n;
      //t = (sphi - delta_y_n / u) * ilphi;
      //y11 = source_y - t;
      y11 = (-sphi + delta_y_n / u) * ilphi;
    }
    int h00 = int(std::floor(y00 * ipix_step + hpix0));
    int h01 = int(std::floor(y01 * ipix_step + hpix0));
    int hmin = std::min(h00, h01);
    int hmax = std::max(h00, h01);
    int h10 = int(std::floor(y10 * ipix_step + hpix0));
    hmin = std::min(hmin, h10);
    hmax = std::max(hmax, h10);
    int h11 = int(std::floor(y11 * ipix_step + hpix0));
    hmin = std::max(std::min(hmin, h11), 0);
    hmax = std::min(std::max(hmax, h11), n_h - 1);
    // If it intercepts voxel then calc alpha_xy_0, alpha_xy_1 and store
    pixel_ptr h_offset = ah_offset + sl_int(hmin) * sl_int(n_v);
    int h_pos = ah_pos + hmin * n_v;
    for (int h = hmin; h <= hmax; h++) {
      const real p2_x = cdetx[a] - sphi * h_pixels[h];
      const real p2_y = sdetx[a] + cphi * h_pixels[h];
      const real delta_x = p2_x - p1_x;
      const real delta_y = p2_y - p1_y;
      if (std::abs(delta_x) < epsilon) {
	report_error("Ooops - delta x");
      } else if (std::abs(delta_y) < epsilon) {
	report_error("Ooops - delta y");
      } else {
	const real inv_dx = 1.0 / delta_x;
	const real inv_dy = 1.0 / delta_y;
	const real alpha_x_0 = delta_x_0 * inv_dx;
	const real alpha_y_0 = delta_y_0 * inv_dy;
	const real alpha_x_n = delta_x_n * inv_dx;
	const real alpha_y_n = delta_y_n * inv_dy;
	const real alpha_x_min = std::min(alpha_x_0, alpha_x_n);
	const real alpha_x_max = std::max(alpha_x_0, alpha_x_n);
	const real alpha_y_min = std::min(alpha_y_0, alpha_y_n);
	const real alpha_y_max = std::max(alpha_y_0, alpha_y_n);
	const real alpha_min = std::max(std::max(alpha_x_min, alpha_y_min),
					real(0.0));
	const real alpha_max = std::min(std::min(alpha_x_max, alpha_y_max),
					real(1.0));
	if (alpha_min < alpha_max - epsilon) {
	  alpha_xy_0[count] = alpha_min;
	  alpha_xy_1[count] = alpha_max;
	  ah_arr[count] = h_offset;
	  ah_index[count] = h_pos;
	  count++;
	}
      }
      h_offset += n_v;
      h_pos += n_v;
    }
    ah_offset += nah;
    ah_pos += nah;
  }
}

void CCPi::cone_beam::b2D_cpu(const real source_x, const real source_y,
			      const real source_z, const real detector_x,
			      const real_1d &h_pixels, const real_1d &v_pixels,
			      const real_1d &angles, pixel_data &pixels,
			      voxel_data &voxels, const int n_angles,
			      const int n_h, const int n_v,
			      const real vox_origin[3], const real vox_size[3],
			      const int nx, const int ny, const int nz,
			      const recon_2d &d_conv)
{
  // Todo - check that source_z to det_z max z angle < 45 for safe usage.
  int mid = -1;
  for (int i = 0; i < n_v; i++) {
    if (v_pixels[i] > source_z) {
      mid = i;
      break;
    }
  }

#pragma omp parallel shared(h_pixels, v_pixels, pixels, voxels, angles, vox_size, vox_origin) firstprivate(n_angles, n_h, n_v, nx, ny, nz)
  {
    int thread_id = omp_get_thread_num();
    int nthreads = omp_get_num_threads();
    real_1d c_angle(n_angles);
    real_1d s_angle(n_angles);
    real_1d p1x(n_angles);
    real_1d p1y(n_angles);
    //real_1d cpy(n_angles);
    //real_1d spy(n_angles);
    real_1d cdetx(n_angles);
    real_1d sdetx(n_angles);
    real_1d ilcphi(n_angles);
    real_1d ilsphi(n_angles);
    const real l = detector_x - source_x;
    for (int a = 0; a < n_angles; a++) {
      real cos_phi = std::cos(angles[a]);
      real sin_phi = std::sin(angles[a]);
      c_angle[a] = cos_phi;
      s_angle[a] = sin_phi;
      real cpy = cos_phi * source_y;
      real spy = sin_phi * source_y;
      p1x[a] = cos_phi * source_x - spy;
      p1y[a] = sin_phi * source_x + cpy;
      cdetx[a] = cos_phi * detector_x;
      sdetx[a] = sin_phi * detector_x;
      ilcphi[a] = l / cos_phi;
      ilsphi[a] = l / sin_phi;
    }

    //const recon_type inv_dz = recon_type(real(1.0) / vox_size[2]);
    const recon_type pzbz = recon_type((source_z - vox_origin[2]) / vox_size[2]);
    const real v_step = v_pixels[1] - v_pixels[0];
    const recon_type pzdv = recon_type((source_z - v_pixels[0]) / v_step);
    const recon_type z_1 = recon_type((vox_origin[2] + vox_size[2]
				       - source_z) / v_step);
    const recon_type z_nm = recon_type((vox_origin[2] + real(nz - 1) * vox_size[2]
					- source_z) / v_step);

    recon_1d delta_z(n_v);
    for (int i = 0; i < n_v; i++)
      delta_z[i] = v_pixels[i] - source_z;
    recon_1d inv_delz(n_v);
    for (int i = 0; i < n_v; i++)
      inv_delz[i] = 1.0 / delta_z[i];
    for (int i = 0; i < n_v; i++)
      delta_z[i] /= vox_size[2];
    recon_1d vox_z(nz + 1);
    for (int i = 0; i <= nz; i++)
      vox_z[i] = vox_origin[2] + real(i) * vox_size[2] - source_z;

    real_1d yvals(ny + 1);
    for (int j = 0; j <= ny; j++)
      yvals[j] = vox_origin[1] + real(j) * vox_size[1];

    int pix_per_vox = 1;
    while (nz * pix_per_vox < n_v)
      pix_per_vox++;
    //const int pix_per_vox = n_v / (nz - 1);
    // How big should the array be
    int xy_size = 2 * pix_per_vox * (n_angles + 10);
    if (xy_size % 32 != 0)
      xy_size += 32 - (xy_size % 32);
    pixel_ptr_1d ah_arr(xy_size);
    int_1d ah_index(xy_size);
    //int_1d h_arr(2 * pix_per_vox * n_angles);
    recon_1d alpha_xy_0(xy_size);
    recon_1d alpha_xy_1(xy_size);
    int nwork = 0;

    const int x_block = nx;
    const int y_block = ny;
    const int a_block = n_angles;
    long counter = 0;
    for (int block_x = 0; block_x < nx; block_x += x_block) {
      int x_step = x_block;
      if (block_x + x_step > nx)
	x_step = nx - block_x;
      for (int block_y = 0; block_y < ny; block_y += y_block) {
	int y_step = y_block;
	if (block_y + y_step > ny)
	  y_step = ny - block_y;
	for (int block_a = 0; block_a < n_angles; block_a += a_block) {
	  int a_step = a_block;
	  if (block_a + a_step > n_angles)
	    a_step = n_angles - block_a;
	  // we don't block h since we have a min/max from the x/y blocks

	  for (int ix = 0; ix < x_step; ix++) {
	    int_1d kv(n_v);
	    int i = block_x + ix;
	    const real x_0 = vox_origin[0] + real(i) * vox_size[0];
	    const real x_n = vox_origin[0] + real(i + 1) * vox_size[0];
	    for (int jx = 0; jx < y_step; jx++) {
	      if (counter % nthreads == thread_id) {
		int j = block_y + jx;
		bproject_ah(source_x, source_y, pixels, voxels, x_0, yvals[j],
			    x_n, yvals[j + 1], nz, i, j, a_step, n_h, n_v,
			    h_pixels, mid, c_angle, s_angle, delta_z, inv_delz,
			    vox_z, pzbz, pzdv, z_1, z_nm, p1x, p1y, cdetx,
			    sdetx, ilcphi, ilsphi, block_a, ah_arr,
			    alpha_xy_0, alpha_xy_1, ah_index, block_a, nwork);
		if (nwork > 0) {
		  if (nwork > xy_size)
		    report_error("back project overflow");
		  else {
#ifdef CBZCHECK
		    // find safe region where all m stay inside the main block
		    // like min_v_all/max_v_all in calc_xy
		    recon_type pzdeps = pzdv * (1.0 + epsilon);
		    int min_v = 0;
		    for (int m = 0; m < count; m++) {
		      int v = int(std::ceil(pzdeps + z_1 / alpha_xy_0[m]));
		      if (v > min_v)
			min_v = v;
		    }
		    if (min_v > 0)
		      std::cerr << "Lower z range wrong for cone-beam "
				<< min_v << '\n';
		    int max_v = n_v;
		    for (int m = 0; m < count; m++) {
		      int v = int(std::floor(pzdv + z_nm / alpha_xy_0[m]));
		      if (v < max_v)
			max_v = v;
		    }
		    if (max_v < n_v)
		      std::cerr << "Upper z range wrong for cone-beam "
				<< max_v << '\n';
#endif // CBZCHECK
		    calc_ah_z(ah_arr, &(voxels[i][j][0]), alpha_xy_0,
			      alpha_xy_1, nwork, pzbz, n_v, nz, mid, delta_z,
			      inv_delz, vox_z, kv);
		  }
		}
	      }
	      counter++;
	    }
	  }
	}
      }
    }      
  }
}

void CCPi::cone_beam::b2D(const real source_x, const real source_y,
			  const real source_z, const real detector_x,
			  const real_1d &h_pixels, const real_1d &v_pixels,
			  const real_1d &angles, pixel_data &pixels,
			  voxel_data &voxels, const int n_angles,
			  const int n_h, const int n_v,
			  const real vox_origin[3], const real vox_size[3],
			  const int nx, const int ny, const int nz,
			  const bool limited_memory)
{
  // path length from source to detector is independent of rotation
  recon_2d d_conv(boost::extents[n_h][n_v]);
  if (limited_memory) {
    recon_2d id_conv(boost::extents[n_h][n_v]);
    real x2 = (detector_x - source_x) * (detector_x - source_x);
#pragma omp parallel for shared(h_pixels, v_pixels) firstprivate(x2, n_h, n_v, source_y, source_z) schedule(dynamic)
    for (int i = 0; i < n_h; i++) {
      real xy2 = x2 + (h_pixels[i] - source_y) * (h_pixels[i] - source_y);
      for (int j = 0; j < n_v; j++) {
	real x3 = std::sqrt(xy2 + (v_pixels[j] - source_z)
			    * (v_pixels[j] - source_z));
	d_conv[i][j] = recon_type(x3);
	id_conv[i][j] = real(1.0) / x3;
      }
    }
#pragma omp parallel for shared(pixels, d_conv), firstprivate(n_angles, n_h, n_v) schedule(dynamic)
    for (int a = 0; a < n_angles; a++) {
      const recon_type *const dc = assume_aligned(&(d_conv[0][0]), recon_type);
      pixel_type *const pix = assume_aligned(&(pixels[a][0][0]), pixel_type);
      for (int hv = 0; hv < n_h * n_v; hv++)
	pix[hv] *= dc[hv];
    }
    if (machine::has_accelerator())
      b2D_accel(source_x, source_y, source_z, detector_x, h_pixels, v_pixels,
		angles, pixels, voxels, n_angles, n_h, n_v, vox_origin,
		vox_size, nx, ny, nz, d_conv);
    else
      b2D_cpu(source_x, source_y, source_z, detector_x, h_pixels, v_pixels,
	      angles, pixels, voxels, n_angles, n_h, n_v, vox_origin,
	      vox_size, nx, ny, nz, d_conv);          
    // This will have some numerical noise.
#pragma omp parallel for shared(pixels, id_conv), firstprivate(n_angles, n_h, n_v) schedule(dynamic)
    for (int a = 0; a < n_angles; a++) {
      const recon_type *const dc = assume_aligned(&(id_conv[0][0]), recon_type);
      pixel_type *const pix = assume_aligned(&(pixels[a][0][0]), pixel_type);
      for (int hv = 0; hv < n_h * n_v; hv++)
	pix[hv] *= dc[hv];
    }
  } else {
    real x2 = (detector_x - source_x) * (detector_x - source_x);
#pragma omp parallel for shared(h_pixels, v_pixels) firstprivate(x2, n_h, n_v, source_y, source_z) schedule(dynamic)
    for (int i = 0; i < n_h; i++) {
      real xy2 = x2 + (h_pixels[i] - source_y) * (h_pixels[i] - source_y);
      for (int j = 0; j < n_v; j++) {
	d_conv[i][j] = std::sqrt(xy2 + (v_pixels[j] - source_z)
				 * (v_pixels[j] - source_z));
      }
    }
    pixel_3d dpixels(boost::extents[n_angles][n_h][n_v]);
#pragma omp parallel for shared(dpixels, pixels, d_conv), firstprivate(n_angles, n_h, n_v) schedule(dynamic)
    for (int a = 0; a < n_angles; a++) {
      const recon_type *const dc = assume_aligned(&(d_conv[0][0]), recon_type);
      const pixel_type *const pix = assume_aligned(&(pixels[a][0][0]),
						   pixel_type);
      pixel_type *const dpx = assume_aligned(&(dpixels[a][0][0]), pixel_type);
      for (int hv = 0; hv < n_h * n_v; hv++)
	dpx[hv] = pix[hv] * dc[hv];
    }
    if (machine::has_accelerator())
      b2D_accel(source_x, source_y, source_z, detector_x, h_pixels, v_pixels,
		angles, dpixels, voxels, n_angles, n_h, n_v, vox_origin,
		vox_size, nx, ny, nz, d_conv);
    else
      b2D_cpu(source_x, source_y, source_z, detector_x, h_pixels, v_pixels,
	      angles, dpixels, voxels, n_angles, n_h, n_v, vox_origin, vox_size,
	      nx, ny, nz, d_conv);
  }
}

#ifdef USE_OPENCL

void CCPi::cone_beam::f2D_accel(const real source_x, const real source_y,
				const real source_z, const real detector_x,
				const real_1d &h_pixels,
				const real_1d &v_pixels, const real_1d &angles,
				pixel_data &pixels, voxel_data &voxels,
				const int n_angles, const int n_h,
				const int n_v, const real grid_offset[3],
				const real voxel_size[3], const int nx_voxels,
				const int ny_voxels, const int nz_voxels)
{
  // Todo - check that source_z to det_z max z angle < 45 for safe usage.
  int mid = -1;
  for (int i = 0; i < n_v; i++) {
    if (v_pixels[i] > source_z) {
      mid = i;
      break;
    }
  }
  if (v_pixels[mid - 1] > source_z - epsilon and
      v_pixels[mid - 1] < source_z + epsilon)
    report_error("Oops");
  if (v_pixels[mid] > source_z - epsilon and
      v_pixels[mid] < source_z + epsilon)
    report_error("Oops");
  int nthreads = machine::number_of_accelerators();

#pragma omp parallel shared(h_pixels, v_pixels, pixels, voxels, angles, voxel_size, grid_offset) firstprivate(n_angles, n_h, n_v, nx_voxels, ny_voxels, nz_voxels, nthreads) num_threads(nthreads)
  {
    int thread_id = omp_get_thread_num();
    const char kernel_name[] = "cone_xy_z";
    if (!machine::check_accelerator_kernel(kernel_name, thread_id))
      report_error("forward kernel problem");
    else {
      // path length from source to detector is independent of rotation
      recon_2d d_conv(boost::extents[n_h][n_v]);
      real x2 = (detector_x - source_x) * (detector_x - source_x);
      for (int i = 0; i < n_h; i++) {
	real xy2 = x2 + (h_pixels[i] - source_y) * (h_pixels[i] - source_y);
	for (int j = 0; j < n_v; j++) {
	  d_conv[i][j] = std::sqrt(xy2 + (v_pixels[j] - source_z)
				   * (v_pixels[j] - source_z));
	}
      }

      const recon_type inv_dz = recon_type(real(1.0) / voxel_size[2]);
      const recon_type pzbz = recon_type((source_z
					  - grid_offset[2]) / voxel_size[2]);
      recon_1d delta_z(n_v);
      for (int i = 0; i < n_v; i++)
	delta_z[i] = v_pixels[i] - source_z;
      recon_1d inv_delz(n_v);
      for (int i = 0; i < n_v; i++)
	inv_delz[i] = 1.0 / delta_z[i];
      recon_1d vox_z(nz_voxels + 32);
      for (int i = 0; i <= nz_voxels; i++)
	vox_z[i] = grid_offset[2] + real(i) * voxel_size[2] - source_z;

      dev_ptr dp_del_z = machine::device_allocate(sl_int(n_v)
						  * sizeof(float),
						  true, thread_id);
      dev_ptr dp_inv_z = machine::device_allocate(sl_int(n_v)
						  * sizeof(float),
						  true, thread_id);
      dev_ptr dp_vox_z = machine::device_allocate(sl_int(nz_voxels + 32)
						  * sizeof(float),
						  true, thread_id);
      machine::copy_to_device(delta_z.data(), dp_del_z,
			      n_v * sizeof(float), thread_id);
      machine::copy_to_device(inv_delz.data(), dp_inv_z,
			      n_v * sizeof(float), thread_id);
      machine::copy_to_device(vox_z.data(), dp_vox_z,
			      (nz_voxels + 32) * sizeof(float), thread_id);

      const real pixel_step = (h_pixels[1] - h_pixels[0]);
      const real ipix_step = 1.0 / pixel_step;
      const real hpix0 = (source_y - h_pixels[0]) / pixel_step;
      const real l = detector_x - source_x;

      // calc space on accelerator for copy
      sl_int accel_size =
	machine::largest_alloc(thread_id) / sizeof(voxel_type);
      double xy_accel = double(accel_size) / double(nz_voxels);
      int proj_x = nx_voxels;
      int proj_y = ny_voxels;
      bool flip = true;
      while (sl_int(proj_x) * sl_int(proj_y) > xy_accel) {
	if (flip)
	  proj_x /= 2;
	else
	  proj_y /= 2;
	flip = !flip;
      }
      const int x_block = proj_x;
      const int y_block = proj_y;

      sl_int nyz = sl_int(ny_voxels) * sl_int(nz_voxels);
      int max_n = std::max(x_block, y_block);
      int ah_size = 2 * max_n + 2;
      if (ah_size % 32 != 0)
	ah_size += 32 - (ah_size % 32);
      recon_1d l_xy(ah_size);
      voxel_ptr_1d ij_arr(ah_size);
      int_1d ij_offsets(ah_size);
      int nwork = 0;

      sl_int accel_proj = n_angles;
      int ah_max = std::max(n_v, ah_size);
      sl_int mem_size = (accel_proj * sl_int(n_h))
	* (sl_int(n_v) * sizeof(pixel_type)
	   + sl_int(ah_size) * (sizeof(recon_type) + sizeof(int)))
	+ sl_int(x_block) * sl_int(y_block) * sl_int(nz_voxels)
	* sizeof(voxel_type);
      // this is really just 4
      sl_int data_max = std::max(std::max(sizeof(voxel_type),
					  sizeof(recon_type)), sizeof(int));
      int factor = 1;
      while (mem_size > machine::available_mem(thread_id) or
	     sl_int(accel_proj) * sl_int(n_h) * sl_int(ah_max)
	     * data_max > machine::largest_alloc(thread_id)) {
	factor++;
	accel_proj = n_angles / factor;
	if (n_angles % factor != 0)
	  accel_proj++;
	mem_size = (accel_proj * sl_int(n_h))
	  * (sl_int(n_v) * sizeof(pixel_type)
	     + sl_int(ah_size) * (sizeof(recon_type) + sizeof(int)))
	  + sl_int(x_block) * sl_int(y_block) * sl_int(nz_voxels)
	  * sizeof(voxel_type);
      }
      const int a_block = accel_proj;

      // Need 3D as can't overwrite 1D ones while doing async copy.
      // n_h or 2 * max(x_block, y_block)?
      boost::multi_array<int, 2> ah3_offsets(boost::extents[a_block * n_h][ah_size]);
      boost::multi_array<recon_type, 2> l3_xy(boost::extents[a_block * n_h][ah_size]);
      boost::multi_array<int, 1> work_sizes(boost::extents[a_block * n_h]);
      boost::multi_array<int, 1> h_arr(boost::extents[a_block * n_h]);
      
      dev_ptr pix_buf = machine::device_allocate(a_block * sl_int(n_h)
						 * sl_int(n_v)
						 * sizeof(pixel_type), false,
						 thread_id);
      dev_ptr vox_buf = machine::device_allocate(sl_int(proj_x)
						 * sl_int(proj_y)
						 * sl_int(nz_voxels)
						 * sizeof(voxel_type), true,
						 thread_id);
      dev_ptr xy_offsets = machine::device_allocate(sl_int(a_block)
						    * sl_int(n_h)
						    * ah_size * sizeof(int),
						    true, thread_id);
      dev_ptr ij_buff = machine::device_allocate(sl_int(a_block)
						 * sl_int(n_h) * ah_size
						 * sizeof(recon_type),
						 true, thread_id);
      dev_ptr ij_work = machine::device_allocate(sl_int(a_block)
						 * sl_int(n_h)
						 * sizeof(int),
						 true, thread_id);
      dev_ptr h_work = machine::device_allocate(sl_int(a_block)
						* sl_int(n_h)
						* sizeof(int),
						true, thread_id);

      long counter = 0;
      event_t vox_ev;
      // we don't block h since we have a min/max from the x/y blocks
      for (int block_x = 0; block_x < nx_voxels; block_x += x_block) {
	int x_step = x_block;
	if (block_x + x_step > nx_voxels)
	  x_step = nx_voxels - block_x;
	sl_int i_base = sl_int(block_x) * nyz;
	real vx = grid_offset[0] + real(block_x) * voxel_size[0];
	real wx = grid_offset[0] + real(block_x + x_step) * voxel_size[0];
	for (int block_y = 0; block_y < ny_voxels; block_y += y_block) {
	  int y_step = y_block;
	  if (block_y + y_step > ny_voxels)
	    y_step = ny_voxels - block_y;
	  machine::copy_to_device(&voxels[block_x][block_y][0], vox_buf,
				  y_step, nz_voxels * sizeof(voxel_type),
				  x_step, y_step, nz_voxels * sizeof(float),
				  thread_id, &vox_ev);
	  sl_int block_yz = y_step * nz_voxels;
	  sl_int x_base = sl_int(block_x) * block_yz;
	  sl_int ij_base = i_base + sl_int(block_y) * sl_int(nz_voxels);
	  real vy = grid_offset[1] + real(block_y) * voxel_size[1];
	  real wy = grid_offset[1] + real(block_y + y_step) * voxel_size[1];
	  for (int block_a = 0; block_a < n_angles; block_a += a_block) {
	    int a_step = a_block;
	    if (block_a + a_step > n_angles)
	      a_step = n_angles - block_a;
	    if (counter % nthreads == thread_id) {
	      event_t ah_ev;
	      machine::copy_to_device(&pixels[block_a][0][0], pix_buf,
				      a_step * sl_int(n_h)
				      * sl_int(n_v)
				      * sizeof(pixel_type), thread_id,
				      &ah_ev);
	      sl_int ah_base = x_base + sl_int(block_y) * sl_int(nz_voxels);
	      std::vector<std::vector<event_t > *> ev(a_step);
	      int cx = 0;
	      for (int ax = 0; ax < a_step; ax++) {
		recon_1d alpha_inv(2 * max_n);
		int_1d kv(n_v);
		int a = block_a + ax;
		real cphi = std::cos(angles[a]);
		real sphi = std::sin(angles[a]);
		real p1_x = cphi * source_x - sphi * source_y;
		real p1_y = sphi * source_x + cphi * source_y;

		// from bproj
		//int hmin = 0;
		//int hmax = n_h - 1;
		/**/
		real y00;
		real y01;
		real y10;
		real y11;
		const real delta_x_0 = vx - p1_x;
		const real delta_y_0 = vy - p1_y;
		const real delta_x_n = wx - p1_x;
		const real delta_y_n = wy - p1_y;
		const real cx_0 = cphi * vx - source_x; // Todo - could pass this in
		const real cx_n = cphi * wx - source_x; // cx_0 + cphi * d_x
		const real sy_0 = sphi * vy; // Todo could pass in as array
		const real sy_n = sphi * wy;
		if (std::abs(sphi) > std::abs(cphi)) {
		  const real ilphi = l / sphi;
		  real u = sy_0 + cx_0;
		  y00 = (cphi - delta_x_0 / u) * ilphi;
		  u = sy_n + cx_0;
		  y01 = (cphi - delta_x_0 / u) * ilphi;
		  u = sy_0 + cx_n;
		  y10 = (cphi - delta_x_n / u) * ilphi;
		  u = sy_n + cx_n;
		  y11 = (cphi - delta_x_n / u) * ilphi;
		} else {
		  real ilphi = l / cphi;
		  real u = cx_0 + sy_0;
		  y00 = (-sphi + delta_y_0 / u) * ilphi;
		  u = cx_0 + sy_n;
		  y01 = (-sphi + delta_y_n / u) * ilphi;
		  u = cx_n + sy_0;
		  y10 = (-sphi + delta_y_0 / u) * ilphi;
		  u = cx_n + sy_n;
		  y11 = (-sphi + delta_y_n / u) * ilphi;
		}
		int h00 = int(std::floor(y00 * ipix_step + hpix0));
		int h01 = int(std::floor(y01 * ipix_step + hpix0));
		int hmin = std::min(h00, h01);
		int hmax = std::max(h00, h01);
		int h10 = int(std::floor(y10 * ipix_step + hpix0));
		hmin = std::min(hmin, h10);
		hmax = std::max(hmax, h10);
		int h11 = int(std::floor(y11 * ipix_step + hpix0));
		hmin = std::max(std::min(hmin, h11), 0);
		hmax = std::min(std::max(hmax, h11), n_h - 1);
		/**/
		// end bproj

		int csize = 0;
		int cstart = cx;
		for (int h = hmin; h <= hmax; h++) {
		  // rotate source and detector positions by current angle
		  real p2_x = cphi * detector_x - sphi * h_pixels[h];
		  real p2_y = sphi * detector_x + cphi * h_pixels[h];
		  fproject_xy(p1_x, p1_y, p2_x, p2_y, pixels, voxels, vx, vy,
			      voxel_size[0], voxel_size[1], x_step, y_step,
			      nz_voxels, a, h, n_v, ij_base, nyz, l_xy, ij_arr,
			      ij_offsets, block_yz, ah_base, nwork);
		  if (nwork > 0) {
		    if (nwork > ah_size) {
		      report_error("forward project overflow");
		      nwork = 0;
		    } else {
		      for (int d = 0; d < nwork; d++)
			ah3_offsets[cx][d] = ij_offsets[d];
		      for (int d = 0; d < nwork; d++)
			l3_xy[cx][d] = l_xy[d];
		      work_sizes[cx] = nwork;
		      h_arr[cx] = ax * n_h + h;
		      cx++;
		      csize++;
		    }
		  }
		}
		if (csize > 0) {
		  ev[ax] = new std::vector<event_t>(6);
		  machine::copy_to_device(&ah3_offsets[cstart][0],
					  xy_offsets,
					  cstart * ah_size * sizeof(int),
					  csize * ah_size * sizeof(int),
					  thread_id, &((*(ev[ax]))[0]));
		  machine::copy_to_device(&l3_xy[cstart][0], ij_buff,
					  cstart * ah_size
					  * sizeof(recon_type),
					  csize * ah_size
					  * sizeof(recon_type),
					  thread_id, &((*(ev[ax]))[1]));
		  machine::copy_to_device(&work_sizes[cstart], ij_work,
					  cstart * sizeof(int),
					  csize * sizeof(int), thread_id,
					  &((*(ev[ax]))[2]));
		  machine::copy_to_device(&h_arr[cstart], h_work,
					  cstart * sizeof(int),
					  csize * sizeof(int), thread_id,
					  &((*(ev[ax]))[3]));
		  (*(ev[ax]))[4] = ah_ev;
		  (*(ev[ax]))[5] = vox_ev;
		  //int pix_offset = (ax * n_h) * n_v;
		  machine::run_cone_xy(kernel_name, pix_buf, vox_buf,
				       ij_buff, xy_offsets, h_work,
				       ij_work, n_v, nz_voxels, cstart,
				       ah_size, pzbz, mid, dp_del_z, inv_dz,
				       dp_inv_z, dp_vox_z, n_v, csize,
				       thread_id, ev[ax]);
		} else
		  ev[ax] = 0;
	      }
	      machine::accelerator_barrier(thread_id);
	      for (int ix = 0; ix < a_step; ix++)
		if (ev[ix] != 0)
		  delete ev[ix];
	      machine::copy_from_device(pix_buf, &pixels[block_a][0][0],
					a_step * sl_int(n_h)
					* sl_int(n_v)
					* sizeof(pixel_type), thread_id);
	      // Todo - improve this
	      for (int ax = 0; ax < a_step; ax++) {
		int a = block_a + ax;
		const recon_type *const dc = assume_aligned(&(d_conv[0][0]),
							    recon_type);
		pixel_type *const pix = assume_aligned(&(pixels[a][0][0]),
						       pixel_type);
		for (int hv = 0; hv < n_h * n_v; hv++)
		  pix[hv] *= dc[hv];
	      }
	    }
	    counter++;
	  }
	}
      }
      machine::device_free(h_work, thread_id);
      machine::device_free(ij_work, thread_id);
      machine::device_free(ij_buff, thread_id);
      machine::device_free(xy_offsets, thread_id);
      machine::device_free(vox_buf, thread_id);
      machine::device_free(pix_buf, thread_id);
      machine::device_free(dp_vox_z, thread_id);
      machine::device_free(dp_inv_z, thread_id);
      machine::device_free(dp_del_z, thread_id);
      machine::accelerator_complete(thread_id);
    }
  }
}

void CCPi::cone_beam::b2D_accel(const real source_x, const real source_y,
				const real source_z, const real detector_x,
				const real_1d &h_pixels,
				const real_1d &v_pixels, const real_1d &angles,
				pixel_data &pixels, voxel_data &voxels,
				const int n_angles, const int n_h,
				const int n_v, const real vox_origin[3],
				const real vox_size[3], const int nx,
				const int ny, const int nz,
				const recon_2d &d_conv)
{
  // Todo - check that source_z to det_z max z angle < 45 for safe usage.
  int mid = -1;
  for (int i = 0; i < n_v; i++) {
    if (v_pixels[i] > source_z) {
      mid = i;
      break;
    }
  }

  int nthreads = machine::number_of_accelerators();

#pragma omp parallel shared(h_pixels, v_pixels, pixels, voxels, angles, vox_size, vox_origin) firstprivate(n_angles, n_h, n_v, nx, ny, nz, nthreads) num_threads(nthreads)
  {
    int thread_id = omp_get_thread_num();
    const char kernel_name[] = "cone_ah_z";
    if (!machine::check_accelerator_kernel(kernel_name, thread_id))
      report_error("backward kernel problem");
    else {
      real_1d c_angle(n_angles);
      real_1d s_angle(n_angles);
      real_1d p1x(n_angles);
      real_1d p1y(n_angles);
      //real_1d cpy(n_angles);
      //real_1d spy(n_angles);
      real_1d cdetx(n_angles);
      real_1d sdetx(n_angles);
      real_1d ilcphi(n_angles);
      real_1d ilsphi(n_angles);
      const real l = detector_x - source_x;
      for (int a = 0; a < n_angles; a++) {
	real cos_phi = std::cos(angles[a]);
	real sin_phi = std::sin(angles[a]);
	c_angle[a] = cos_phi;
	s_angle[a] = sin_phi;
	real cpy = cos_phi * source_y;
	real spy = sin_phi * source_y;
	p1x[a] = cos_phi * source_x - spy;
	p1y[a] = sin_phi * source_x + cpy;
	cdetx[a] = cos_phi * detector_x;
	sdetx[a] = sin_phi * detector_x;
	ilcphi[a] = l / cos_phi;
	ilsphi[a] = l / sin_phi;
      }

      const recon_type inv_dz = recon_type(real(1.0) / vox_size[2]);
      const recon_type pzbz = recon_type((source_z
					  - vox_origin[2]) / vox_size[2]);
      const real v_step = v_pixels[1] - v_pixels[0];
      const recon_type pzdv = recon_type((source_z - v_pixels[0]) / v_step);
      const recon_type z_1 = recon_type((vox_origin[2] + vox_size[2]
					 - source_z) / v_step);
      const recon_type z_nm = recon_type((vox_origin[2] + real(nz - 1)
					  * vox_size[2] - source_z) / v_step);

      recon_1d delta_z(n_v);
      for (int i = 0; i < n_v; i++)
	delta_z[i] = v_pixels[i] - source_z;
      recon_1d inv_delz(n_v);
      for (int i = 0; i < n_v; i++)
	inv_delz[i] = 1.0 / delta_z[i];
      for (int i = 0; i < n_v; i++)
	delta_z[i] *= inv_dz;
      recon_1d vox_z(nz + 32);
      for (int i = 0; i <= nz; i++)
	vox_z[i] = vox_origin[2] + real(i) * vox_size[2] - source_z;

      dev_ptr dp_del_z = machine::device_allocate(sl_int(n_v)
						  * sizeof(float),
						  true, thread_id);
      dev_ptr dp_inv_z = machine::device_allocate(sl_int(n_v)
						  * sizeof(float),
						  true, thread_id);
      dev_ptr dp_vox_z = machine::device_allocate(sl_int(nz + 32)
						  * sizeof(float),
						  true, thread_id);
      machine::copy_to_device(delta_z.data(), dp_del_z,
			      n_v * sizeof(float), thread_id);
      machine::copy_to_device(inv_delz.data(), dp_inv_z,
			      n_v * sizeof(float), thread_id);
      machine::copy_to_device(vox_z.data(), dp_vox_z,
			      (nz + 32) * sizeof(float), thread_id);

      real_1d yvals(ny + 1);
      for (int j = 0; j <= ny; j++)
	yvals[j] = vox_origin[1] + real(j) * vox_size[1];

      // calc space on accelerator for copy
      sl_int accel_size =
	machine::largest_alloc(thread_id) / sizeof(pixel_type);
      int accel_proj = n_angles;
      int factor = 1;
      while (accel_proj * sl_int(n_h) * sl_int(n_v) >accel_size) {
	factor++;
	accel_proj = n_angles / factor;
	if (n_angles % factor != 0)
	  accel_proj++;
      }
      const int a_block = accel_proj;

      int pix_per_vox = 1;
      while (nz * pix_per_vox < n_v)
	pix_per_vox++;
      //const int pix_per_vox = n_v / (nz - 1);
      // How big should the array be
      int xy_size = 2 * pix_per_vox * (n_angles + 10);
      if (xy_size % 32 != 0)
	xy_size += 32 - (xy_size % 32);
      pixel_ptr_1d ah_arr(xy_size);
      int_1d ah_index(xy_size);
      //int_1d h_arr(2 * pix_per_vox * n_angles);
      recon_1d alpha_xy_0(xy_size);
      recon_1d alpha_xy_1(xy_size);
      int nwork = 0;

      int x_block = nx;
      int y_block = ny;
      int xy_max = std::max(nz, xy_size);
      sl_int mem_size = accel_proj * sl_int(n_h) * sl_int(n_v)
	* sizeof(pixel_type)
	+ sl_int(x_block) * sl_int(y_block) * sl_int(xy_max)
	* (sizeof(voxel_type) + sizeof(recon_type) + sizeof(int));
      // this is really just 4
      sl_int data_max = std::max(std::max(sizeof(voxel_type),
					  sizeof(recon_type)), sizeof(int));
      bool flip = true;
      while (mem_size > machine::available_mem(thread_id) or
	     sl_int(x_block) * sl_int(y_block) * sl_int(xy_max) * data_max >
	     machine::largest_alloc(thread_id)) {
	if (flip)
	  x_block /= 2;
	else
	  y_block /= 2;
	flip = !flip;
	mem_size = accel_proj * sl_int(n_h) * sl_int(n_v)
	  * sizeof(pixel_type)
	  + sl_int(x_block) * sl_int(y_block) * sl_int(xy_max)
	  * (sizeof(voxel_type) + 2 * sizeof(recon_type) + sizeof(int));
      }
      /*
	add_output("b2D sizes ");
	add_output(a_block);
	add_output(" ");
	add_output(x_block);
	add_output(" ");
	add_output(y_block);
	add_output(" ");
	add_output(xy_size);
	send_output();
      */

      // Need 3D as can't overwrite 1D ones while doing async copy.
      boost::multi_array<int, 2> ah3_offsets(boost::extents[x_block * y_block][xy_size]);
      boost::multi_array<recon_type, 2> l3_xy_0(boost::extents[x_block * y_block][xy_size]);
      boost::multi_array<recon_type, 2> l3_xy_1(boost::extents[x_block * y_block][xy_size]);
      boost::multi_array<int, 1> work_sizes(boost::extents[x_block * y_block]);
      boost::multi_array<int, 1> h_arr(boost::extents[x_block * y_block]);
      
      dev_ptr pix_buf = machine::device_allocate(accel_proj
						 * sl_int(n_h)
						 * sl_int(n_v)
						 * sizeof(pixel_type), true,
						 thread_id);
      dev_ptr vox_buf = machine::device_allocate(sl_int(x_block)
						 * sl_int(y_block)
						 * sl_int(nz)
						 * sizeof(voxel_type), false,
						 thread_id);
      dev_ptr xy_offsets = machine::device_allocate(sl_int(x_block)
						    * sl_int(y_block)
						    * xy_size * sizeof(int),
						    true, thread_id);
      dev_ptr xy0_buff = machine::device_allocate(sl_int(x_block)
						  * sl_int(y_block) * xy_size
						  * sizeof(recon_type),
						  true, thread_id);
      dev_ptr xy1_buff = machine::device_allocate(sl_int(x_block)
						  * sl_int(y_block) * xy_size
						  * sizeof(recon_type),
						  true, thread_id);
      dev_ptr xy_work = machine::device_allocate(sl_int(x_block)
						 * sl_int(y_block)
						 * sizeof(int),
						 true, thread_id);
      dev_ptr h_work = machine::device_allocate(sl_int(x_block)
						* sl_int(y_block)
						* sizeof(int),
						true, thread_id);

      int_1d kv(n_v);
      long counter = 0;
      for (int ap = 0; ap < n_angles; ap += accel_proj) {
	int p_step = accel_proj;
	if (ap + p_step > n_angles)
	  p_step = n_angles - ap;
	event_t pixel_ev;
	machine::copy_to_device(&pixels[ap][0][0], pix_buf,
				p_step * sl_int(n_h) * sl_int(n_v)
				* sizeof(pixel_type), thread_id, &pixel_ev);
	for (int block_x = 0; block_x < nx; block_x += x_block) {
	  int x_step = x_block;
	  if (block_x + x_step > nx)
	    x_step = nx - block_x;
	  for (int block_y = 0; block_y < ny; block_y += y_block) {
	    int y_step = y_block;
	    if (block_y + y_step > ny)
	      y_step = ny - block_y;
	    if (counter % nthreads == thread_id) {
	      event_t vox_ev;
	      machine::copy_to_device(&voxels[block_x][block_y][0], vox_buf,
				      ny, nz * sizeof(voxel_type),
				      x_step, y_step, nz * sizeof(float),
				      thread_id, &vox_ev);
	      for (int block_a = 0; block_a < n_angles; block_a += a_block) {
		int a_step = a_block;
		if (block_a + a_step > n_angles)
		  a_step = n_angles - block_a;
		int offset = 0;
		std::vector<std::vector<event_t> *> ev(x_block);
		// we don't block h since we have a min/max from the x/y blocks

		int csize = 0;
		int cstart = 0;
		for (int ix = 0; ix < x_step; ix++) {
		  int i = block_x + ix;
		  const real x_0 = vox_origin[0] + real(i) * vox_size[0];
		  const real x_n = vox_origin[0] + real(i + 1) * vox_size[0];
		  for (int jx = 0; jx < y_step; jx++) {
		    int j = block_y + jx;
		    bproject_ah(source_x, source_y, pixels, voxels, x_0,
				yvals[j], x_n, yvals[j + 1], nz, i, j, a_step,
				n_h, n_v, h_pixels, mid, c_angle, s_angle,
				delta_z, inv_delz, vox_z, pzbz, pzdv, z_1, z_nm,
				p1x, p1y, cdetx, sdetx, ilcphi, ilsphi,
				ap + block_a, ah_arr, alpha_xy_0, alpha_xy_1,
				ah_index, block_a, nwork);
		    offset += xy_size;
		    if (nwork > 0) {
		      if (nwork > xy_size) {
			report_error("back project overflow");
			nwork = 0;
		      } else {
			for (int d = 0; d < nwork; d++)
			  ah3_offsets[csize][d] = ah_index[d];
			for (int d = 0; d < nwork; d++)
			  l3_xy_0[csize][d] = alpha_xy_0[d];
			for (int d = 0; d < nwork; d++)
			  l3_xy_1[csize][d] = alpha_xy_1[d];
			h_arr[csize] = ix * y_step + jx;
			work_sizes[csize] = nwork;
			csize++;
		      }
		    }
		  }
		  if (csize > 0) {
		    ev[ix] = new std::vector<event_t>(7);
		    machine::copy_to_device(&ah3_offsets[cstart][0], xy_offsets,
					    cstart * xy_size * sizeof(int),
					    csize * xy_size * sizeof(int),
					    thread_id, &((*(ev[ix]))[0]));
		    machine::copy_to_device(&l3_xy_0[cstart][0], xy0_buff,
					    cstart * xy_size
					    * sizeof(recon_type),
					    csize * xy_size
					    * sizeof(recon_type),
					    thread_id, &((*(ev[ix]))[1]));
		    machine::copy_to_device(&l3_xy_1[cstart][0], xy1_buff,
					    cstart * xy_size
					    * sizeof(recon_type),
					    csize * xy_size
					    * sizeof(recon_type),
					    thread_id, &((*(ev[ix]))[2]));
		    machine::copy_to_device(&work_sizes[cstart], xy_work,
					    cstart * sizeof(int),
					    csize * sizeof(int), thread_id,
					    &((*(ev[ix]))[3]));
		    machine::copy_to_device(&h_arr[cstart], h_work,
					    cstart * sizeof(int),
					    csize * sizeof(int), thread_id,
					    &((*(ev[ix]))[4]));
		    (*(ev[ix]))[5] = pixel_ev;
		    (*(ev[ix]))[6] = vox_ev;
		    //int vox_offset = (ix * y_step) * nz;
		    machine::run_cone_ah(kernel_name, pix_buf, vox_buf,
					 xy0_buff, xy1_buff, xy_offsets,
					 xy_work, h_work, n_v, nz, cstart,
					 xy_size, pzbz, mid, dp_del_z, inv_dz,
					 dp_inv_z, dp_vox_z, n_v, csize,
					 thread_id, ev[ix]);
		  } else
		    ev[ix] = 0;
		}
		machine::accelerator_barrier(thread_id);
		for (int ix = 0; ix < x_step; ix++)
		  if (ev[ix] != 0)
		    delete ev[ix];
	      }
	      // no events copy is blocking
	      machine::copy_from_device(&voxels[block_x][block_y][0], vox_buf,
					ny, nz * sizeof(voxel_type),
					x_step, y_step, nz * sizeof(float),
					thread_id);
	    }
	    counter++;
	  }
	}
      }
    }
  }      
}

#else

void CCPi::cone_beam::f2D_accel(const real source_x, const real source_y,
				const real source_z, const real detector_x,
				const real_1d &h_pixels,
				const real_1d &v_pixels, const real_1d &angles,
				pixel_data &pixels, voxel_data &voxels,
				const int n_angles, const int n_h,
				const int n_v, const real grid_offset[3],
				const real voxel_size[3], const int nx_voxels,
				const int ny_voxels, const int nz_voxels)
{
  report_error("BUG: attempt to run on accelerator");
}

void CCPi::cone_beam::b2D_accel(const real source_x, const real source_y,
				const real source_z, const real detector_x,
				const real_1d &h_pixels,
				const real_1d &v_pixels, const real_1d &angles,
				pixel_data &pixels, voxel_data &voxels,
				const int n_angles, const int n_h,
				const int n_v, const real vox_origin[3],
				const real vox_size[3], const int nx,
				const int ny, const int nz,
				const recon_2d &d_conv)
{
  report_error("BUG: attempt to run on accelerator");
}

#endif // OPENCL
