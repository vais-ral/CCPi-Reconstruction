
#include <map>
#include <vector>
#include <cmath>
#include <cstring>
#include <float.h>
#ifdef MATLAB_MEX_FILE
#  include "mex_types.hpp"
#else
#  include "base_types.hpp"
#endif // mex
#include "instruments.hpp"
#include "ui_calls.hpp"
#include "accel.hpp"
#ifdef TEST2D
#  include <iostream>
#endif // TEST2D
#if defined(__AVX__)
#  include <immintrin.h>
#endif // AVX etc
#include <omp.h>

static const recon_type epsilon = FLT_EPSILON;

void CCPi::parallel_beam::gen_mapping(int_1d &mapping, int &map_type,
				      const real_1d &v_pixels, const real vox_z,
				      const real size_z, const int nv)
{
  for (int v = 0; v < nv; v++)
    mapping[v] = std::floor((v_pixels[v] - vox_z) / size_z);
  map_type = 0;
  bool found = true;
  for (int i = 0; i < nv; i++) {
    if (mapping[i] != i) {
      found = false;
      break;
    }
  }
  // test for 1,2,4 pixels_per_voxel, else just use mapping array k
  if (found)
    map_type = 1;
  else {
    found = true;
    for (int i = 0; i < nv; i++) {
      if (mapping[i] != i / 2) {
	found = false;
	break;
      }
    }
    if (found)
      map_type = 2;
    else {
      found = true;
      for (int i = 0; i < nv; i++) {
	if (mapping[i] != i / 4) {
	  found = false;
	  break;
	}
      }
      if (found)
	map_type = 4;
    }
  }
}

void CCPi::parallel_beam::calc_xy_z(pixel_type *const pixels,
				    const voxel_ptr_1d &voxels,
				    const recon_1d &l_xy, const int n,
				    const int nv, const int nz,
				    const int_1d &mapping, const int map_type)
{
  pixel_type *const pix = assume_aligned(pixels, pixel_type);
  recon_type *lptr = assume_aligned(&(l_xy[0]), recon_type);
  switch (map_type) {
  case 1:
#if defined(__AVX2__) && PIXEL_SIZE == 4 && !defined(__GNUC__)
    for (int m = 0; m < n; m++) {
      const voxel_type *const vox = voxels[m];
      const recon_type alpha = lptr[m];
      __m256 al = _mm256_set1_ps(alpha);
      for (int v = 0; v < nv; v += 8) {
	_mm256_store_ps(&pix[v],
			_mm256_fmadd_ps(al,_mm256_load_ps(&vox[v]),
					_mm256_load_ps(&pix[v])));
      }
    }
#elif defined(__AVX__) && PIXEL_SIZE == 4 && !defined(__GNUC__)
    for (int m = 0; m < n; m++) {
      const voxel_type *const vox = voxels[m];
      const recon_type alpha = lptr[m];
      __m256 al = _mm256_set1_ps(alpha);
      for (int v = 0; v < nv; v += 8) {
	_mm256_store_ps(&pix[v],
			_mm256_add_ps(_mm256_mul_ps(al,_mm256_load_ps(&vox[v])),
				      _mm256_load_ps(&pix[v])));
      }
    }
#else
    for (int m = 0; m < n; m++) {
      const voxel_type *const vox = assume_aligned(voxels[m], voxel_type);
      const recon_type alpha = lptr[m];
      for (int v = 0; v < nv; v++)
	pix[v] += vox[v] * alpha;
    }
#endif
    break;
  case 2:
#if defined(__AVX2__) && PIXEL_SIZE == 4
    {
      // assume this is low to high abcdabcd->aabbccdd
      const __m256i offsets = _mm256_set_epi32(0, 4, 1, 5, 2, 6, 3, 7);
      for (int m = 0; m < n; m++) {
	const voxel_type *const vox = voxels[m];
	// __mm256 al = _mm256_broadcast_ss(lptr[m]); ?
	const recon_type alpha = lptr[m];
	__m256 al = _mm256_set1_ps(alpha);
	int h = 0;
	for (int l = 0; l < nv; l += 8) {
	  __m256 vx = _mm256_broadcast_ps((__m128 *)&vox[h]); //abcdabcd
	  __m256 v = _mm256_permutevar8x32_ps(vx, offsets);
	  _mm256_store_ps(&pix[l], _mm256_fmadd_ps(al, v,
						   _mm256_load_ps(&pix[l])));
	  h += 4;
	}
      }
    }
#elif defined(__AVX__) && PIXEL_SIZE == 4
    for (int m = 0; m < n; m++) {
      const voxel_type *const vox = voxels[m];
      // __mm256 al = _mm256_broadcast_ss(lptr[m]); ?
      const recon_type alpha = lptr[m];
      __m256 al = _mm256_set1_ps(alpha);
      int h = 0;
      // Todo - can load 8 way vox as 2x means 16 way pix unroll
      for (int l = 0; l < nv; l += 8) {
	__m256 vx = _mm256_broadcast_ps((__m128 *)&vox[h]); //abcdabcd
	__m256 vl = _mm256_unpacklo_ps(vx, vx);
	__m256 vh = _mm256_unpackhi_ps(vx, vx);
	__m256 v = _mm256_blend_ps(vl, vh, 0xf0); //aabbccdd I hope
	_mm256_store_ps(&pix[l], _mm256_add_ps(_mm256_mul_ps(al, v),
					       _mm256_load_ps(&pix[l])));
	h += 4;
      }
    }
#else
    for (int m = 0; m < n; m++) {
      const voxel_type *const vox = assume_aligned(voxels[m], voxel_type);
      const recon_type alpha = lptr[m];
      int v = 0;
      for (int l = 0; l < nz; l++) {
	pix[v + 0] += vox[l] * alpha;
	pix[v + 1] += vox[l] * alpha;
	v += 2;
      }
    }
#endif
    break;
  case 4:
#if defined(__AVX2__) && PIXEL_SIZE == 4
    for (int m = 0; m < n; m++) {
      const voxel_type *const vox = voxels[m];
      const recon_type alpha = lptr[m];
      __m256 al = _mm256_set1_ps(alpha);
      // Todo - try to load a vector and rearrange?
      int h = 0;
      for (int l = 0; l < nv; l += 8) {
	__m256 vl = _mm256_broadcast_ss(&vox[h + 0]); //aaaaaaaa
	__m256 vh = _mm256_broadcast_ss(&vox[h + 1]); //bbbbbbbb
	__m256 v = _mm256_blend_ps(vl, vh, 0xf0); //aaaabbbb
	_mm256_store_ps(&pix[l], _mm256_fmadd_ps(al, v,
						 _mm256_load_ps(&pix[l])));
	h += 2;
      }
    }
#elif defined(__AVX__) && PIXEL_SIZE == 4
    for (int m = 0; m < n; m++) {
      const voxel_type *const vox = voxels[m];
      const recon_type alpha = lptr[m];
      __m256 al = _mm256_set1_ps(alpha);
      // Todo - try to load a vector and rearrange?
      int h = 0;
      for (int l = 0; l < nv; l += 8) {
	__m256 vl = _mm256_broadcast_ss(&vox[h + 0]); //aaaaaaaa
	__m256 vh = _mm256_broadcast_ss(&vox[h + 1]); //bbbbbbbb
	__m256 v = _mm256_blend_ps(vl, vh, 0xf0); //aaaabbbb
	_mm256_store_ps(&pix[l], _mm256_add_ps(_mm256_mul_ps(al, v),
					       _mm256_load_ps(&pix[l])));
	h += 2;
      }
    }
#else
    for (int m = 0; m < n; m++) {
      const voxel_type *const vox = assume_aligned(voxels[m], voxel_type);
      const recon_type alpha = lptr[m];
      int v = 0;
      for (int l = 0; l < nz; l++) {
	pix[v + 0] += vox[l] * alpha;
	pix[v + 1] += vox[l] * alpha;
	pix[v + 2] += vox[l] * alpha;
	pix[v + 3] += vox[l] * alpha;
	v += 4;
      }
    }
#endif
    break;
  default:
    for (int m = 0; m < n; m++) {
      const voxel_type *const vox = assume_aligned(voxels[m], voxel_type);
      const recon_type alpha = lptr[m];
      for (int v = 0; v < nv; v++)
	pix[v] += vox[mapping[v]] * alpha;
    }
    break;
  }
}

void CCPi::parallel_beam::fproject_xy(const real p2_x, const real p2_y,
				      pixel_data &pixels, voxel_data &voxels,
				      const real b_x, const real b_y,
				      const real d_x, const real d_y,
				      const int nx, const int ny, const int nz,
				      const int a, const int h, const int nv,
				      const recon_type d_conv, const real cphi,
				      const real sphi, const sl_int ij_base,
				      const sl_int nyz, recon_1d &l_xy,
				      voxel_ptr_1d &ij_arr, int_1d &ij_index,
				      const sl_int block_yz, const int xy_base,
				      int &count)
{
  // Todo? In parallel beam some of this should be common in a since
  // all h within a have same angle to voxels.
  count = 0;
  //sl_int nyz = sl_int(ny) * sl_int(nz);
  if (std::abs(cphi) < epsilon) {
    if (std::abs(sphi) < epsilon) {
      // Its not a line - shouldn't happen
      //return;
    } else {
      // line parallel to y
      int i = int(std::floor((p2_x - b_x) / d_x)); // == q2_x
      if (i >= 0 and i < nx) {
	if (sphi < 0.0) {
	  voxel_ptr ij_offset = voxels.data() + ij_base + sl_int(i) * nyz
	    + sl_int(ny - 1) * sl_int(nz);
	  int xy_offset = xy_base + i * block_yz + (ny - 1) * nz;
	  for (int j = ny - 1; j >= 0; j--) {
	    l_xy[count] = d_y;
	    ij_arr[count] = ij_offset;
	    ij_index[count] = xy_offset;
	    ij_offset -= nz;
	    xy_offset -= nz;
	    count++;
	  }
	} else {
	  voxel_ptr ij_offset = voxels.data() + ij_base + sl_int(i) * nyz;
	  int xy_offset = xy_base + i * block_yz;
	  for (int j = 0; j < ny; j++) {
	    l_xy[count] = d_y;
	    ij_arr[count] = ij_offset;
	    ij_index[count] = xy_offset;
	    ij_offset += nz;
	    xy_offset += nz;
	    count++;
	  }
	}
      }
    }
  } else if (std::abs(sphi) < epsilon) {
    // line parallel to x
    int j = int(std::floor((p2_y - b_y) / d_y)); // == q2_y
    if (j >= 0 and j < ny) {
      if (cphi < 0.0) {
	voxel_ptr ij_offset = voxels.data() + ij_base + sl_int(nx - 1) * nyz
	  + sl_int(j) * sl_int(nz);
	int xy_offset = xy_base + (nx - 1) * block_yz + j * nz;
	for (int i = nx - 1; i >= 0; i--) {
	  l_xy[count] = d_x;
	  ij_arr[count] = ij_offset;
	  ij_index[count] = xy_offset;
	  ij_offset -= nyz;
	  xy_offset -= block_yz;
	  count++;
	}
      } else {
	voxel_ptr ij_offset = voxels.data() + ij_base + sl_int(j) * sl_int(nz);
	int xy_offset = xy_base + j * nz;
	for (int i = 0; i < nx; i++) {
	  l_xy[count] = d_x;
	  ij_arr[count] = ij_offset;
	  ij_index[count] = xy_offset;
	  ij_offset += nyz;
	  xy_offset += block_yz;
	  count++;
	}
      }
    }
  } else {
    // general line drawing - find intercepts
    const real x_0 = b_x;
    const real y_0 = b_y;
    const real x_n = b_x + real(nx) * d_x;
    const real y_n = b_y + real(ny) * d_y;
    const real delta_x = cphi;
    const real delta_y = sphi;
    const real inv_dx = 1.0 / delta_x;
    const real inv_dy = 1.0 / delta_y;
    const real alpha_x_0 = (x_0 - p2_x) * inv_dx;
    const real alpha_y_0 = (y_0 - p2_y) * inv_dy;
    const real alpha_x_n = (x_n - p2_x) * inv_dx;
    const real alpha_y_n = (y_n - p2_y) * inv_dy;
    const real alpha_x_min = std::min(alpha_x_0, alpha_x_n);
    const real alpha_x_max = std::max(alpha_x_0, alpha_x_n);
    const real alpha_y_min = std::min(alpha_y_0, alpha_y_n);
    const real alpha_y_max = std::max(alpha_y_0, alpha_y_n);
    const real alpha_min = std::max(std::max(alpha_x_min, alpha_y_min),
				    real(-d_conv));
    const real alpha_max = std::min(std::min(alpha_x_max, alpha_y_max),
				    real(0.0));
    if (alpha_min < alpha_max - epsilon) {
      real_1d alpha_x(nx + 1);
      for (int i = 0; i <= nx; i++)
	alpha_x[i] = (b_x + real(i) * d_x - p2_x) * inv_dx;
      real_1d alpha_y(ny + 1);
      for (int i = 0; i <= ny; i++)
	alpha_y[i] = (b_y + real(i) * d_y - p2_y) * inv_dy;
      int x = 0;
      int y = 0;
      real alpha_p = alpha_min;
      if (delta_x > 0.0) {
	if (delta_y > 0.0) {
	  if (alpha_min == alpha_x_0) {
	    x = 0;
	    y = int(std::floor((p2_y + alpha_min * delta_y - b_y) / d_y));
	  } else if (alpha_min == alpha_y_0) {
	    x = int(std::floor((p2_x + alpha_min * delta_x - b_x) / d_x));
	    y = 0;
	  } else
	    report_error("something wrong in x+ y+");
	  // could do x_next/y_next here and only calc the one that changes
	  // inside the if statements below, which would reduce the flops
	  voxel_ptr ij_offset = voxels.data() + ij_base + sl_int(x) * nyz
	    + sl_int(y) * sl_int(nz);
	  int xy_offset = xy_base + x * block_yz + y * nz;
	  while (x < nx and y < ny) {
	    ij_arr[count] = ij_offset;
	    ij_index[count] = xy_offset;
	    if (alpha_x[x + 1] < alpha_y[y + 1] - epsilon) {
	      l_xy[count] = (alpha_x[x + 1] - alpha_p);
	      alpha_p = alpha_x[x + 1];
	      ij_offset += nyz;
	      xy_offset += block_yz;
	      x++;
	    } else if (alpha_y[y + 1] < alpha_x[x + 1] - epsilon) {
	      l_xy[count] = (alpha_y[y + 1] - alpha_p);
	      alpha_p = alpha_y[y + 1];
	      ij_offset += nz;
	      xy_offset += nz;
	      y++;
	    } else {
	      real mx = std::max(alpha_x[x + 1], alpha_y[y + 1]);
	      l_xy[count] = (mx - alpha_p);
	      x++;
	      y++;
	      ij_offset += nyz + nz;
	      xy_offset += block_yz + nz;
	      alpha_p = mx;
	    }
	    count++;
	  }
	} else {
	  if (alpha_min == alpha_y_n) {
	    x = int(std::floor((p2_x + alpha_min * delta_x - b_x) / d_x));
	    y = ny - 1;
	  } else if (alpha_min == alpha_x_0) {
	    x = 0;
	    y = int(std::floor((p2_y + alpha_min * delta_y - b_y) / d_y));
	  } else
	    report_error("something wrong in x+ y-");
	  voxel_ptr ij_offset = voxels.data() + ij_base + sl_int(x) * nyz
	    + sl_int(y) * sl_int(nz);
	  int xy_offset = xy_base + x * block_yz + y * nz;
	  while (x < nx and y >= 0) {
	    ij_arr[count] = ij_offset;
	    ij_index[count] = xy_offset;
	    if (alpha_x[x + 1] < alpha_y[y] - epsilon) {
	      l_xy[count] = (alpha_x[x + 1] - alpha_p);
	      alpha_p = alpha_x[x + 1];
	      ij_offset += nyz;
	      xy_offset += block_yz;
	      x++;
	    } else if (alpha_y[y] < alpha_x[x + 1] - epsilon) {
	      l_xy[count] = (alpha_y[y] - alpha_p);
	      alpha_p = alpha_y[y];
	      ij_offset -= nz;
	      xy_offset -= nz;
	      y--;
	    } else {
	      real mx = std::max(alpha_x[x + 1], alpha_y[y]);
	      l_xy[count] = (mx - alpha_p);
	      x++;
	      y--;
	      ij_offset += nyz - nz;
	      xy_offset += block_yz - nz;
	      alpha_p = mx;
	    }
	    count++;
	  }
	}
      } else {
	if (delta_y > 0.0) {
	  if (alpha_min == alpha_x_n) {
	    x = nx - 1;
	    y = int(std::floor((p2_y + alpha_min * delta_y - b_y) / d_y));
	  } else if (alpha_min == alpha_y_0) {
	    x = int(std::floor((p2_x + alpha_min * delta_x - b_x) / d_x));
	    y = 0;
	  } else
	    report_error("something wrong in x- y+");
	  voxel_ptr ij_offset = voxels.data() + ij_base + sl_int(x) * nyz
	    + sl_int(y) * sl_int(nz);
	  int xy_offset = xy_base + x * block_yz + y * nz;
	  while (x >= 0 and y < ny) {
	    ij_arr[count] = ij_offset;
	    ij_index[count] = xy_offset;
	    if (alpha_x[x] < alpha_y[y + 1] - epsilon) {
	      l_xy[count] = (alpha_x[x] - alpha_p);
	      alpha_p = alpha_x[x];
	      ij_offset -= nyz;
	      xy_offset -= block_yz;
	      x--;
	    } else if (alpha_y[y + 1] < alpha_x[x] - epsilon) {
	      l_xy[count] = (alpha_y[y + 1] - alpha_p);
	      alpha_p = alpha_y[y + 1];
	      ij_offset += nz;
	      xy_offset += nz;
	      y++;
	    } else {
	      real mx = std::max(alpha_y[y + 1], alpha_x[x]);
	      l_xy[count] = (mx - alpha_p);
	      x--;
	      y++;
	      ij_offset += -nyz + nz;
	      xy_offset += -block_yz + nz;
	      alpha_p = mx;
	    }
	    count++;
	  }
	} else {
	  if (alpha_min == alpha_x_n) {
	    x = nx - 1;
	    if (alpha_min == alpha_y_n)
	      y = ny - 1;
	    else
	      y = int(std::floor((p2_y + alpha_min * delta_y - b_y) / d_y));
	  } else if (alpha_min == alpha_y_n) {
	    x = int(std::floor((p2_x + alpha_min * delta_x - b_x) / d_x));
	    y = ny - 1;
	  } else
	    report_error("something wrong in x- y-");
	  voxel_ptr ij_offset = voxels.data() + ij_base + sl_int(x) * nyz
	    + sl_int(y) * sl_int(nz);
	  int xy_offset = xy_base + x * block_yz + y * nz;
	  while (x >= 0 and y >= 0) {
	    ij_arr[count] = ij_offset;
	    ij_index[count] = xy_offset;
	    if (alpha_x[x] < alpha_y[y] - epsilon) {
	      l_xy[count] = (alpha_x[x] - alpha_p);
	      alpha_p = alpha_x[x];
	      ij_offset -= nyz;
	      xy_offset -= block_yz;
	      x--;
	    } else if (alpha_y[y] < alpha_x[x] - epsilon) {
	      l_xy[count] = (alpha_y[y] - alpha_p);
	      alpha_p = alpha_y[y];
	      ij_offset -= nz;
	      xy_offset -= nz;
	      y--;
	    } else {
	      real mx = std::max(alpha_x[x], alpha_y[y]);
	      l_xy[count] = (mx - alpha_p);
	      x--;
	      y--;
	      ij_offset -= (nyz + nz);
	      xy_offset -= (block_yz + nz);
	      alpha_p = mx;
	    }
	    count++;
	  }
	}
      }
    }
  }
}

void CCPi::parallel_beam::f2D_cpu(const real_1d &h_pixels,
				  const real_1d &v_pixels,
				  const real_1d &angles, const int n_angles,
				  const int nh_pixels, const int nv_pixels,
				  const real vox_origin[3],
				  const real vox_size[3],
				  const int nx, const int ny, const int nz,
				  pixel_data &pixels, voxel_data &voxels)
{
  // set detector z to 2* the yz limits of the voxels, so it misses
  // longest voxel dim should be sqrt(3), so 2 should be safe
  real detector_x = real(2.0) * std::max(std::abs(vox_origin[0]),
					 std::max(std::abs(vox_origin[1]),
						  std::abs(vox_origin[2])));
  // at 0 degrees
  //real p2_x = detector_x;
  //real p2_y = 0.0;
  //real p1_x = p2_x - real(3.0) * detector_x = -2.0 * detector_x
  //real p1_y = 0.0;
  // path length from source to detector is independent of rotation
  recon_type d_conv = recon_type(3.0 * detector_x);

  int_1d mapping(nv_pixels);
  int map_type = 0;
  gen_mapping(mapping, map_type, v_pixels, vox_origin[2], vox_size[2],
	      nv_pixels);

  const real ihp_step = 1.0 / (h_pixels[1] - h_pixels[0]);
  const real h_pix0 = h_pixels[0] / (h_pixels[1] - h_pixels[0]);
  sl_int nyz = sl_int(ny) * sl_int(nz);
  int max_n = std::max(nx, ny);
  recon_1d l_xy(2 * max_n);
  voxel_ptr_1d ij_arr(2 * max_n + 1);
  int_1d ij_offsets(2 * max_n + 1);
  int nwork = 0;

  const int a_block = n_angles;
  const int x_block = 32;
  const int y_block = 32;
  for (int block_a = 0; block_a < n_angles; block_a += a_block) {
    int a_step = a_block;
    if (block_a + a_step > n_angles)
      a_step = n_angles - block_a;
    // we don't block h since we have a min/max from the x/y blocks
    for (int block_x = 0; block_x < nx; block_x += x_block) {
      int x_step = x_block;
      if (block_x + x_step > nx)
	x_step = nx - block_x;
      sl_int i_base = sl_int(block_x) * nyz;
      real vx = vox_origin[0] + real(block_x) * vox_size[0];
      real wx = vox_origin[0] + real(block_x + x_step) * vox_size[0];
      for (int block_y = 0; block_y < ny; block_y += y_block) {
	int y_step = y_block;
	if (block_y + y_step > ny)
	  y_step = ny - block_y;
	sl_int ij_base = i_base + sl_int(block_y) * sl_int(nz);
	real vy = vox_origin[1] + real(block_y) * vox_size[1];
	real wy = vox_origin[1] + real(block_y + y_step) * vox_size[1];

#pragma omp parallel for shared(h_pixels, v_pixels, pixels, voxels, angles, d_conv, vox_size, vox_origin, mapping) firstprivate(n_angles, nh_pixels, nv_pixels, nx, ny, nz, detector_x, ihp_step, h_pix0, nyz, ij_base, vx, vy, wx, wy, map_type) schedule(dynamic)
	for (int ax = 0; ax < a_step; ax++) {
	  int a = block_a + ax;
	  // rotate source and detector positions by current angle
	  real cphi = std::cos(angles[a]);
	  real sphi = std::sin(angles[a]);

          // from bproject
          int hmin = 0;
          int hmax = nh_pixels - 1;
          if (std::abs(cphi) < epsilon) {
            // line is parallel to y, so just find h that is between x_0/x_n
            if (sphi < 0.0) {
              hmin = int(std::floor(vx * ihp_step - h_pix0));
              if (hmin < 0)
                hmin = 0;
              // decrease x_n by tol so inside upper boundary, not on it?
              hmax = int(std::floor((wx - epsilon) * ihp_step - h_pix0));
              if (hmax >= nh_pixels)
                hmax = nh_pixels - 1;
            } else {
              hmin = int(std::ceil((- wx + epsilon) * ihp_step - h_pix0));
              if (hmin < 0)
                hmin = 0;
              // decrease x_n by tol so inside upper boundary, not on it?
              hmax = int(std::ceil((- vx) * ihp_step - h_pix0));
              if (hmax >= nh_pixels)
                hmax = nh_pixels - 1;
            }
          } else if (std::abs(sphi) < epsilon) {
            if (cphi < 0.0) {
              hmin = int(std::ceil((- wy + epsilon) * ihp_step - h_pix0));
              if (hmin < 0)
                hmin = 0;
              // decrease wy by tol so inside upper boundary, not on it?
              hmax = int(std::ceil((- vy) * ihp_step - h_pix0));
              if (hmax >= nh_pixels)
                hmax = nh_pixels - 1;
            } else {
              hmin = int(std::floor(vy * ihp_step - h_pix0));
              if (hmin < 0)
                hmin = 0;
              // decrease wy by tol so inside upper boundary, not on it?
              hmax = int(std::floor((wy - epsilon) * ihp_step - h_pix0));
              if (hmax >= nh_pixels)
                hmax = nh_pixels - 1;
             }
          } else {
            real ymin;
            real ymax;
            // Todo - pass in arrays of vy * cphi etc?
            if (cphi > 0.0) {
              if (sphi > 0.0) {
                ymin = vy * cphi - wx * sphi; //y10;
                ymax = wy * cphi - vx * sphi; //y01;
              } else {
                ymin = vy * cphi - vx * sphi; //y00;
                ymax = wy * cphi - wx * sphi; //y11;
              }
            } else {
              if (sphi > 0.0) {
                ymin = wy * cphi - wx * sphi; //y11;
                ymax = vy * cphi - vx * sphi; //y00;
              } else {
                ymin = wy * cphi - vx * sphi; //y01;
                ymax = vy * cphi - wx * sphi; //y10;
              }
            }
            if (ymin < h_pixels[0])
              hmin = 0;
            else
              hmin = int(std::floor(ymin * ihp_step - h_pix0));
            //if (ymax >= h_pixels[nh]) Todo
            hmax = std::min(int(std::floor(ymax * ihp_step - h_pix0)),
                            nh_pixels-1);
          }
          // end bproject

	  for (int h = hmin; h <= hmax; h++) {
	    real p2_x = cphi * detector_x - sphi * h_pixels[h];
	    real p2_y = sphi * detector_x + cphi * h_pixels[h];
	    //real p1_x = p2_x - real(3.0) * cphi * detector_x;
	    //real p1_y = p2_y - real(3.0) * sphi * detector_x;
	    fproject_xy(p2_x, p2_y, pixels, voxels, vx, vy, vox_size[0],
			vox_size[1], x_step, y_step, nz, a, h, nv_pixels,
			d_conv, cphi, sphi, ij_base, nyz, l_xy, ij_arr,
			ij_offsets, nz, nz, nwork);
	    if (nwork > 0) {
	      if (nwork > 2 * max_n + 1)
		report_error("forward project overflow");
	      else
		calc_xy_z(&(pixels[a][h][0]), ij_arr, l_xy, nwork, nv_pixels,
			  nz, mapping, map_type);
	    }
	  }
	}
      }
    }
  }
}

void CCPi::parallel_beam::f2D(const real_1d &h_pixels, const real_1d &v_pixels,
			      const real_1d &angles, const int n_angles,
			      const int nh_pixels, const int nv_pixels,
			      const real vox_origin[3], const real vox_size[3],
			      const int nx, const int ny, const int nz,
			      pixel_data &pixels, voxel_data &voxels)
{
  if (machine::has_accelerator())
    f2D_accel(h_pixels, v_pixels, angles, n_angles, nh_pixels, nv_pixels,
	      vox_origin, vox_size, nx, ny, nz, pixels, voxels);
  else
    f2D_cpu(h_pixels, v_pixels, angles, n_angles, nh_pixels, nv_pixels,
	    vox_origin, vox_size, nx, ny, nz, pixels, voxels);
}

void CCPi::parallel_beam::calc_ah_z(const pixel_ptr_1d &pixels,
				    voxel_type *const voxels,
				    const recon_1d &l_xy, const int n,
				    const int nv, const int nz,
				    const int_1d &mapping, const int map_type)
{
  voxel_type *const vox = assume_aligned(voxels, voxel_type);
  recon_type *lptr = assume_aligned(&(l_xy[0]), recon_type);
  switch (map_type) {
  case 1:
#if defined(__AVX2__) && PIXEL_SIZE == 4 && !defined(__GNUC__)
    for (int m = 0; m < n; m++) {
      const pixel_type *const pix = pixels[m];
      const recon_type alpha = lptr[m];
      __m256 al = _mm256_set1_ps(alpha);
      for (int v = 0; v < nv; v += 8) {
	_mm256_store_ps(&vox[v],
			_mm256_fmadd_ps(al,_mm256_load_ps(&pix[v]),
					_mm256_load_ps(&vox[v])));
      }
    }
#elif defined(__AVX__) && PIXEL_SIZE == 4 && !defined(__GNUC__)
    for (int m = 0; m < n; m++) {
      const pixel_type *const pix = pixels[m];
      const recon_type alpha = lptr[m];
      __m256 al = _mm256_set1_ps(alpha);
      for (int v = 0; v < nv; v += 8) {
	_mm256_store_ps(&vox[v],
			_mm256_add_ps(_mm256_mul_ps(al,_mm256_load_ps(&pix[v])),
				      _mm256_load_ps(&vox[v])));
      }
    }
#else
    for (int m = 0; m < n; m++) {
      const pixel_type *const pix = assume_aligned(pixels[m], pixel_type);
      const recon_type alpha = lptr[m];
      for (int v = 0; v < nv; v++)
	vox[v] += pix[v] * alpha;
    }
#endif
    break;
  case 2:
#if defined(__AVX2__) && PIXEL_SIZE == 4
    {
      //const __m256i offsets = _mm256_set_epi32(0, 2, 4, 6, 1, 3, 6, 7);
      // Todo - use 8 way vox vectors
      for (int m = 0; m < n; m++) {
	const pixel_type *const pix = pixels[m];
	const recon_type alpha = lptr[m];
	__m128 al = _mm_set1_ps(alpha);
	int v = 0;
	for (int l = 0; l < nz; l += 4) {
	  __m256 p = _mm256_load_ps(&pix[v]);
	  __m128 pl = _mm256_extractf128_ps(p, 0x0); //abcd
	  __m128 ph = _mm256_extractf128_ps(p, 0x1); //efgh
	  _mm_store_ps(&vox[l], _mm_fmadd_ps(al, _mm_hadd_ps(pl, ph),
					     _mm_load_ps(&vox[l])));
	  /*
	  __m128 pp = _mm256_permutevar8x32_ps(p, offsets); //acegbdfh
	  __m128 p0 = _mm256_extractf128_ps(pp, 0x0); //aceg
	  __m128 p1 = _mm256_extractf128_ps(pp, 0x1); //bdfh
	  // sum pair of values since both a multiplied by alpha
	  _mm_store_ps(&vox[l], _mm_fmadd_ps(al, _mm_add_ps(p0, p1),
	           _mm_load_ps(&vox[l])));
	  */
	  v += 8;
	}
      }
    }
#elif defined(__AVX__) && PIXEL_SIZE == 4
    for (int m = 0; m < n; m++) {
      const pixel_type *const pix = pixels[m];
      const recon_type alpha = lptr[m];
      __m256 al = _mm256_set1_ps(alpha);
      int v = 0;
      for (int l = 0; l < nz; l += 8) {
	// opt manual recommends reducing port-5 pressure with insert
	__m256 pl = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_load_ps(&pix[v + 0])),
					 _mm_load_ps(&pix[v + 8]), 0x1);
	__m256 ph = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_load_ps(&pix[v + 4])),
					 _mm_load_ps(&pix[v + 12]), 0x1);
	/*
	__m256 p0 = _mm256_load_ps(&pix[v + 0]); //abcdefgh
	__m256 p1 = _mm256_load_ps(&pix[v + 8]); //ijklmnop
	__m256 pl = _mm256_permute2f128_ps(p0, p1, 0x20); //abcdijkl
	__m256 ph = _mm256_permute2f128_ps(p0, p1, 0x31); //efghmnop
	*/
	// hadd should now give (ab)(cd)(ef)(gh)(ij)(kl)(mn)(op)
	_mm256_store_ps(&vox[l],
			_mm256_add_ps(_mm256_mul_ps(al, _mm256_hadd_ps(pl, ph)),
				      _mm256_load_ps(&vox[l])));
	v += 16;
      }
    }
#else
    for (int m = 0; m < n; m++) {
      const pixel_type *const pix = assume_aligned(pixels[m], pixel_type);
      const recon_type alpha = lptr[m];
      int v = 0;
      for (int l = 0; l < nz; l++) {
	vox[l] += (pix[v + 0] + pix[v + 1]) * alpha;
	v += 2;
      }
    }
#endif
    break;
  case 4:
#if defined(__AVX2__) && PIXEL_SIZE == 4
    for (int m = 0; m < n; m++) {
      const pixel_type *const pix = pixels[m];
      const recon_type alpha = lptr[m];
      __m128 al = _mm_set1_ps(alpha);
      int v = 0;
      for (int l = 0; l < nz; l += 4) {
	__m256 p0 = _mm256_load_ps(&pix[v + 0]);
	__m256 p1 = _mm256_load_ps(&pix[v + 8]);
	__m256 p = _mm256_hadd_ps(p0, p1);
	__m128 pl = _mm256_extractf128_ps(p, 0x0); //abcd
	__m128 ph = _mm256_extractf128_ps(p, 0x1); //efgh
	_mm_store_ps(&vox[l], _mm_fmadd_ps(al, _mm_hadd_ps(pl, ph),
					   _mm_load_ps(&vox[l])));
	v += 16;
      }
    }
#elif defined(__AVX__) && PIXEL_SIZE == 4
    for (int m = 0; m < n; m++) {
      const pixel_type *const pix = pixels[m];
      const recon_type alpha = lptr[m];
      __m256 al = _mm256_set1_ps(alpha);
      int v = 0;
      for (int l = 0; l < nz; l += 8) {
	__m256 pa = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_load_ps(&pix[v + 0])),
					 _mm_load_ps(&pix[v + 16]), 0x1);
	__m256 pb = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_load_ps(&pix[v + 4])),
					 _mm_load_ps(&pix[v + 20]), 0x1);
	__m256 pc = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_load_ps(&pix[v + 8])),
					 _mm_load_ps(&pix[v + 24]), 0x1);
	__m256 pd = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_load_ps(&pix[v + 12])),
					 _mm_load_ps(&pix[v + 28]), 0x1);
	/*
	__m256 p0 = _mm256_load_ps(&pix[v + 0]);
	__m256 p1 = _mm256_load_ps(&pix[v + 8]);
	__m256 p2 = _mm256_load_ps(&pix[v + 16]);
	__m256 p3 = _mm256_load_ps(&pix[v + 24]);
	__m256 pa = _mm256_permute2f128_ps(p0, p2, 0x20);
	__m256 pb = _mm256_permute2f128_ps(p0, p2, 0x31);
	__m256 pc = _mm256_permute2f128_ps(p1, p3, 0x20);
	__m256 pd = _mm256_permute2f128_ps(p1, p3, 0x31);
	*/
	__m256 pl = _mm256_hadd_ps(pa, pb);
	__m256 ph = _mm256_hadd_ps(pc, pd);
	_mm256_store_ps(&vox[l],
			_mm256_add_ps(_mm256_mul_ps(al, _mm256_hadd_ps(pl, ph)),
				      _mm256_load_ps(&vox[l])));
	v += 32;
      }
    }
#else
    for (int m = 0; m < n; m++) {
      const pixel_type *const pix = assume_aligned(pixels[m], pixel_type);
      const recon_type alpha = lptr[m];
      int v = 0;
      for (int l = 0; l < nz; l++) {
	vox[l] += (pix[v + 0] + pix[v + 1] + pix[v + 2] + pix[v + 3]) * alpha;
	v += 4;
      }
    }
#endif
    break;
  default:
    for (int m = 0; m < n; m++) {
      const pixel_type *const pix = assume_aligned(pixels[m], pixel_type);
      const recon_type alpha = lptr[m];
      for (int v = 0; v < nv; v++)
	vox[mapping[v]] += pix[v] * alpha;
    }
    break;
  }
}

void CCPi::parallel_beam::bproject_ah(pixel_data &pixels, voxel_data &voxels,
				      const real x_0, const real y_0,
				      const real x_n, const real y_n,
				      const real d_x, const real d_y,
				      const int nz, const int i, const int j,
				      const int n_angles, const int n_h,
				      const int n_v, const real_1d &h_pixels,
				      const real_1d &cangle,
				      const real_1d &sangle,
				      const real_1d &y_offset,
				      const real_1d &i_offset,
				      const real_1d &length, const real h_pix0,
				      const real ihp_step, const int a_off,
				      pixel_ptr_1d &ah_arr, recon_1d &l_xy,
				      int_1d &ah_index, const sl_int a_base,
				      int &count)
{
  // Rather than using the centre just calculate for all 4 corners,
  // generate h values and loop from smallest to largest.
  count = 0;
  // corners (x0,y0), (x0,yn), (xn,y0), (xn,yn)
  // Todo - in parallel we can probably make a better guess at which 2 corners
  // we need for the upper and lower limits.
  //real pixel_step = h_pixels[1] - h_pixels[0];
  sl_int nah = sl_int(n_h) * sl_int(n_v);
  pixel_ptr ah_offset = pixels.data() + sl_int(a_off) * nah;
  int ah_pos = a_base * nah;
  for (int ax = 0; ax < n_angles; ax++) {
    int a = a_off + ax;
    real cphi = cangle[a];
    real sphi = sangle[a];
    // Todo - is there a direction issue here that should swap 0/n?
    if (std::abs(cphi) < epsilon) {
      // line is parallel to y, so just find h that is between x_0/x_n
      if (sphi < 0.0) {
	int hmin = int(std::floor(x_0 * ihp_step - h_pix0));
	if (hmin < 0)
	  hmin = 0;
	// decrease x_n by tol so inside upper boundary, not on it?
	int hmax = int(std::floor((x_n - epsilon) * ihp_step - h_pix0));
	if (hmax >= n_h)
	  hmax = n_h - 1;
	pixel_ptr h_offset = ah_offset + sl_int(hmin) * sl_int(n_v);
	int h_pos = ah_pos + hmin * n_v;
	for (int h = hmin; h <= hmax; h++) {
	  if (h_pixels[h] >= x_0 and h_pixels[h] < x_n) {
	    l_xy[count] = d_y;
	    ah_arr[count] = h_offset;
	    ah_index[count] = h_pos;
	    count++;
	  }
	  h_offset += n_v;
	  h_pos += n_v;
	}
      } else {
	int hmin = int(std::ceil((- x_n + epsilon) * ihp_step - h_pix0));
	if (hmin < 0)
	  hmin = 0;
	// decrease x_n by tol so inside upper boundary, not on it?
	int hmax = int(std::ceil((- x_0) * ihp_step - h_pix0));
	if (hmax >= n_h)
	  hmax = n_h - 1;
	pixel_ptr h_offset = ah_offset + sl_int(hmin) * sl_int(n_v);
	int h_pos = ah_pos + hmin * n_v;
	for (int h = hmin; h <= hmax; h++) {
	  if (-h_pixels[h] >= x_0 and -h_pixels[h] < x_n) {
	    l_xy[count] = d_y;
	    ah_arr[count] = h_offset;
	    ah_index[count] = h_pos;
	    count++;
	  }
	  h_offset += n_v;
	  h_pos += n_v;
	}
      }
    } else if (std::abs(sphi) < epsilon) {
      if (cphi < 0.0) {
	int hmin = int(std::ceil((- y_n + epsilon) * ihp_step - h_pix0));
	if (hmin < 0)
	  hmin = 0;
	// decrease y_n by tol so inside upper boundary, not on it?
	int hmax = int(std::ceil((- y_0) * ihp_step - h_pix0));
	if (hmax >= n_h)
	  hmax = n_h - 1;
	pixel_ptr h_offset = ah_offset + sl_int(hmin) * sl_int(n_v);
	int h_pos = ah_pos + hmin * n_v;
	for (int h = hmin; h <= hmax; h++) {
	  if (- h_pixels[h] >= y_0 and - h_pixels[h] < y_n) {
	    l_xy[count] = d_x;
	    ah_arr[count] = h_offset;
	    ah_index[count] = h_pos;
	    count++;
	  }
	  h_offset += n_v;
	  h_pos += n_v;
	}
      } else {
	int hmin = int(std::floor(y_0 * ihp_step - h_pix0));
	if (hmin < 0)
	  hmin = 0;
	// decrease y_n by tol so inside upper boundary, not on it?
	int hmax = int(std::floor((y_n - epsilon) * ihp_step - h_pix0));
	if (hmax >= n_h)
	  hmax = n_h - 1;
	pixel_ptr h_offset = ah_offset + sl_int(hmin) * sl_int(n_v);
	int h_pos = ah_pos + hmin * n_v;
	for (int h = hmin; h <= hmax; h++) {
	  if (h_pixels[h] >= y_0 and h_pixels[h] < y_n) {
	    l_xy[count] = d_x;
	    ah_arr[count] = h_offset;
	    ah_index[count] = h_pos;
	    count++;
	  }
	  h_offset += n_v;
	  h_pos += n_v;
	}
      }
    } else {
      int hmin;
      int hmax;
      /*
      real y00 = qpx_0 * sphi - qpy_0 * cphi;
      real y01 = qpx_0 * sphi - qpy_n * cphi;
      real y10 = qpx_n * sphi - qpy_0 * cphi;
      real y11 = qpx_n * sphi - qpy_n * cphi;
      */
      real ymin;
      real ymax;
      // Todo - pass in arrays of y_0 * cphi etc?
      if (cphi > 0.0) {
	if (sphi > 0.0) {
	  ymin = y_0 * cphi - x_n * sphi; //y10;
	  ymax = y_n * cphi - x_0 * sphi; //y01;
	} else {
	  ymin = y_0 * cphi - x_0 * sphi; //y00;
	  ymax = y_n * cphi - x_n * sphi; //y11;
	}
      } else {
	if (sphi > 0.0) {
	  ymin = y_n * cphi - x_n * sphi; //y11;
	  ymax = y_0 * cphi - x_0 * sphi; //y00;
	} else {
	  ymin = y_n * cphi - x_0 * sphi; //y01;
	  ymax = y_0 * cphi - x_n * sphi; //y10;
	}
      }
      if (ymin < h_pixels[0])
	hmin = 0;
      else
	hmin = int(std::floor(ymin * ihp_step - h_pix0));
      //if (ymax >= h_pixels[nh]) Todo
      hmax = std::min(int(std::floor(ymax * ihp_step - h_pix0)), n_h-1);
      real ybot = ymin + y_offset[a];
      real ytop = ymax - y_offset[a];
      const real l = length[a];
      pixel_ptr h_offset = ah_offset + sl_int(hmin) * sl_int(n_v);
      int h_pos = ah_pos + hmin * n_v;
      for (int h = hmin; h <= hmax; h++) {
	if (h_pixels[h] > ymin + epsilon) {
	  if (h_pixels[h] < ybot) {
	    l_xy[count] = l * (h_pixels[h] - ymin) * i_offset[a];
	    ah_arr[count] = h_offset;
	    ah_index[count] = h_pos;
	    count++;
	  } else if (h_pixels[h] <= ytop) {
	    l_xy[count] = l;
	    ah_arr[count] = h_offset;
	    ah_index[count] = h_pos;
	    count++;
	  } else if (h_pixels[h] < ymax - epsilon) {
	    l_xy[count] = l * (ymax - h_pixels[h]) * i_offset[a];
	    ah_arr[count] = h_offset;
	    ah_index[count] = h_pos;
	    count++;
	  } else
	    break;
	}
	h_offset += n_v;
	h_pos += n_v;
      }
    }
    ah_offset += nah;
    ah_pos += nah;
  }
}

void CCPi::parallel_beam::b2D_cpu(const real_1d &h_pixels,
				  const real_1d &v_pixels,
				  const real_1d &angles, pixel_data &pixels,
				  voxel_data &voxels, const int n_angles,
				  const int nh_pixels, const int nv_pixels,
				  const real vox_origin[3],
				  const real vox_size[3],
				  const int nx, const int ny, const int nz)
{
  // Todo - 1d arrays of x,y,p1x.p1y positions and 2d p2x/p2y?

  //const real tol = 1e-8;
  //recon_2d line_angles(boost::extents[n_angles][n_h]);
  real_1d c_angle(n_angles);
  real_1d s_angle(n_angles);
  for (int a = 0; a < n_angles; a++) {
    real cos_phi = std::cos(angles[a]);
    real sin_phi = std::sin(angles[a]);
    c_angle[a] = cos_phi;
    s_angle[a] = sin_phi;
  }

  // set detector z to 2* the yz limits of the voxels, so it misses
  // longest voxel dim should be sqrt(3), so 2 should be safe
  //real detector_x = real(2.0) * std::max(std::abs(vox_origin[0]),
  //				 std::max(std::abs(vox_origin[1]),
  //					  std::abs(vox_origin[2])));
  // at 0 degrees
  //real p2_x = detector_x;
  //real p2_y = 0.0;
  // path length from source to detector is independent of rotation
  //recon_type d_conv = recon_type(3.0 * detector_x);
  //const real source_x = -2.0 * detector_x;

  // Todo - do we need this?
  recon_1d vox_z(nz + 1);
  for (int i = 0; i <= nz; i++)
    vox_z[i] = vox_origin[2] + real(i) * vox_size[2];
  real_1d yvals(ny + 1);
  for (int j = 0; j <= ny; j++)
    yvals[j] = vox_origin[1] + real(j) * vox_size[1];

  int_1d mapping(nv_pixels);
  int map_type = 0;
  gen_mapping(mapping, map_type, v_pixels, vox_origin[2], vox_size[2],
	      nv_pixels);

  // rotation of voxel for offsets is same for all voxels
  real_1d y_offset(n_angles);
  real_1d i_offset(n_angles);
  real_1d length(n_angles);
  for (int a = 0; a < n_angles; a++) {
    real cphi = c_angle[a];
    real sphi = s_angle[a];
    // use values for 0,0 voxel x_0 = b_x, x_1 = b_x + d_x etc
    real y00 = vox_origin[1] * cphi - vox_origin[0] * sphi;
    real y01 = (vox_origin[1] + vox_size[1]) * cphi - vox_origin[0] * sphi;
    real y10 = vox_origin[1] * cphi - (vox_origin[0] + vox_size[0]) * sphi;
    real y11 = (vox_origin[1] + vox_size[1]) * cphi
      - (vox_origin[0] + vox_size[0]) * sphi;
    real ymin;
    real ybot;
    if (cphi > 0.0) {
      if (sphi > 0.0) {
	ymin = y10;
	if (y11 < y00) {
	  ybot = y11;
	} else {
	  ybot = y00;
	}
      } else {
	ymin = y00;
	if (y01 < y10) {
	  ybot = y01;
	} else {
	  ybot = y10;
	}
      }
    } else {
      if (sphi > 0.0) {
	ymin = y11;
	if (y01 < y10) {
	  ybot = y01;
	} else {
	  ybot = y10;
	}
      } else {
	ymin = y01;
	if (y11 < y00) {
	  ybot = y11;
	} else {
	  ybot = y00;
	}
      }
    }
    y_offset[a] = ybot - ymin;
    i_offset[a] = 1.0 / y_offset[a];
    if (std::abs(cphi) > std::abs(sphi))
      length[a] = vox_size[0] / std::abs(cphi);
    else
      length[a] = vox_size[1] / std::abs(sphi);
  }
  const real ihp_step = 1.0 / (h_pixels[1] - h_pixels[0]);
  const real h_pix0 = h_pixels[0] / (h_pixels[1] - h_pixels[0]);
  // from bproject_ah
  const int pix_per_vox = nv_pixels / (nz - 1);
  // How big should the array be - Todo - use mapping for pix_per_vox?
  int nwork = 0;
  pixel_ptr_1d ah_arr(2 * pix_per_vox * (n_angles + 10));
  int_1d ah_offsets(2 * pix_per_vox * (n_angles + 10));
  recon_1d l_xy(2 * pix_per_vox * (n_angles + 10));

  const int x_block = 32;
  const int y_block = 32;
  const int a_block = 40;
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

#pragma omp parallel for shared(h_pixels, v_pixels, pixels, voxels, angles, vox_size, vox_origin, yvals, mapping) firstprivate(n_angles, nh_pixels, nv_pixels, nx, ny, nz, y_offset, i_offset, length, h_pix0, ihp_step, c_angle, s_angle, block_x, block_y, x_step, y_step, a_step, block_a, map_type) schedule(dynamic)
	for (int ix = 0; ix < x_step; ix++) {
	  int i = block_x + ix;
	  const real x_0 = vox_origin[0] + real(i) * vox_size[0];
	  const real x_n = vox_origin[0] + real(i + 1) * vox_size[0];
	  for (int jx = 0; jx < y_step; jx++) {
	    int j = block_y + jx;
	    bproject_ah(pixels, voxels, x_0, yvals[j],
			x_n, yvals[j + 1], vox_size[0], vox_size[1], nz,
			i, j, a_step, nh_pixels, nv_pixels, h_pixels,
			c_angle, s_angle, y_offset, i_offset, length,
			h_pix0, ihp_step, block_a, ah_arr, l_xy, ah_offsets,
			nv_pixels, nwork);
	    if (nwork > 0) {
	      if (nwork > 2 * pix_per_vox * (n_angles + 10))
		report_error("back project overflow");
	      else
		calc_ah_z(ah_arr, &(voxels[i][j][0]), l_xy, nwork,
			  nv_pixels, nz, mapping, map_type);
	    }
	  }
	}
      }
    }
  }
}

void CCPi::parallel_beam::b2D(const real_1d &h_pixels, const real_1d &v_pixels,
			      const real_1d &angles, pixel_data &pixels,
			      voxel_data &voxels, const int n_angles,
			      const int nh_pixels, const int nv_pixels,
			      const real vox_origin[3], const real vox_size[3],
			      const int nx, const int ny, const int nz)
{
  if (machine::has_accelerator())
    b2D_accel(h_pixels, v_pixels, angles, pixels, voxels, n_angles, nh_pixels,
	      nv_pixels, vox_origin, vox_size, nx, ny, nz);
  else
    b2D_cpu(h_pixels, v_pixels, angles, pixels, voxels, n_angles, nh_pixels,
	    nv_pixels, vox_origin, vox_size, nx, ny, nz);
}

#ifdef USE_OPENCL

void CCPi::parallel_beam::f2D_accel(const real_1d &h_pixels,
				    const real_1d &v_pixels,
				    const real_1d &angles, const int n_angles,
				    const int nh_pixels, const int nv_pixels,
				    const real vox_origin[3],
				    const real vox_size[3],
				    const int nx, const int ny, const int nz,
				    pixel_data &pixels, voxel_data &voxels)
{
  int_1d mapping(nv_pixels);
  int map_type = 0;
  gen_mapping(mapping, map_type, v_pixels, vox_origin[2], vox_size[2],
	      nv_pixels);
  if (map_type != 1 and map_type != 2)
    f2D_cpu(h_pixels, v_pixels, angles, n_angles, nh_pixels, nv_pixels,
	    vox_origin, vox_size, nx, ny, nz, pixels, voxels);
  else {
    // run a thread for each device
    int nthreads = machine::number_of_accelerators();
#pragma omp parallel shared(h_pixels, v_pixels, pixels, voxels, angles, vox_size, vox_origin, mapping) firstprivate(n_angles, nh_pixels, nv_pixels, nx, ny, nz,map_type) num_threads(nthreads)
    {
      int thread_id = omp_get_thread_num();
      char kernel_name[32];
      if (map_type == 2)
	strcpy(kernel_name, "parallel_xy_z2");
      else
	strcpy(kernel_name, "parallel_xy_z1");
      if (!machine::check_accelerator_kernel(kernel_name, thread_id))
	report_error("forward kernel problem");
      else {
	real detector_x =
	  real(2.0) * std::max(std::abs(vox_origin[0]),
			       std::max(std::abs(vox_origin[1]),
					std::abs(vox_origin[2])));
	recon_type d_conv = recon_type(3.0 * detector_x);
      
	const real ihp_step = 1.0 / (h_pixels[1] - h_pixels[0]);
	const real h_pix0 = h_pixels[0] / (h_pixels[1] - h_pixels[0]);
	sl_int nyz = sl_int(ny) * sl_int(nz);
	int max_n = std::max(nx, ny);
	int xy_size = 2 * max_n + 2;
	recon_1d l_xy(xy_size);
	voxel_ptr_1d ij_arr(xy_size);
	int_1d ij_offsets(xy_size);
	int nwork = 0;
      
	// calc space on accelerator for copy
	sl_int accel_size =
	  machine::largest_alloc(thread_id) / sizeof(voxel_type);
	double xy_accel = double(accel_size) / double(nz);
	int xy_proj = (int)std::sqrt(xy_accel);
	if (xy_proj > max_n)
	  xy_proj = max_n;

	accel_size = machine::largest_alloc(thread_id) / sizeof(pixel_type);
	sl_int accel_proj = accel_size / (nh_pixels * nv_pixels);
	if (accel_proj > n_angles)
	  accel_proj = n_angles;

	const int a_block = accel_proj;
	const int x_block = 32;
	const int y_block = 32;

	dev_ptr pix_buf = machine::device_allocate(a_block * sl_int(nh_pixels)
						   * sl_int(nv_pixels)
						   * sizeof(pixel_type), false,
						   thread_id);
	dev_ptr vox_buf = machine::device_allocate(sl_int(xy_proj)
						   * sl_int(xy_proj)
						   * sl_int(nz)
						   * sizeof(voxel_type), true,
						   thread_id);
	dev_ptr xy_offsets = machine::device_allocate(xy_size * sizeof(int),
						      true, thread_id);
	dev_ptr xy_buff = machine::device_allocate(xy_size * sizeof(recon_type),
						   true, thread_id);
	long counter = 0;
	for (int xp = 0; xp < nx; xp += xy_proj) {
	  int x_proj = xy_proj;
	  if (xp + x_proj > nx)
	    x_proj = nx - xp;
	  for (int yp = 0; yp < ny; yp += xy_proj) {
	    int y_proj = xy_proj;
	    if (yp + y_proj > ny)
	      y_proj = ny - yp;
	    // Its not contiguous, since not a full y range
	    for (int ix = 0; ix < x_proj; ix++)
	      machine::copy_to_device(&voxels[xp][yp][0],
				      vox_buf, sl_int(ix * y_proj)
				      * sl_int(nz) * sizeof(voxel_type),
				      sl_int(y_proj)
				      * sl_int(nz) * sizeof(voxel_type),
				      thread_id);
	    sl_int block_yz = y_proj * nz;
	    for (int block_a = 0; block_a < n_angles; block_a += a_block) {
	      int a_step = a_block;
	      if (block_a + a_step > n_angles)
		a_step = n_angles - block_a;
	      if (counter % nthreads == thread_id) {
		// Todo - do we really want to copy all?
		machine::copy_to_device(&pixels[block_a][0][0], pix_buf,
					a_step * sl_int(nh_pixels)
					* sl_int(nv_pixels)
					* sizeof(pixel_type), thread_id);
		// we don't block h since we have a min/max from the x/y blocks
		for (int block_x = 0; block_x < x_proj; block_x += x_block) {
		  int x_step = x_block;
		  if (block_x + x_step > x_proj)
		    x_step = x_proj - block_x;
		  sl_int i_base = sl_int(xp + block_x) * nyz;
		  sl_int x_base = sl_int(block_x) * block_yz;
		  real vx = vox_origin[0] + real(xp + block_x) * vox_size[0];
		  real wx = vox_origin[0] + real(xp + block_x + x_step)
		    * vox_size[0];
		  for (int block_y = 0; block_y < y_proj; block_y += y_block) {
		    int y_step = y_block;
		    if (block_y + y_step > y_proj)
		      y_step = y_proj - block_y;
		    sl_int ij_base = i_base + sl_int(yp + block_y) * sl_int(nz);
		    real vy = vox_origin[1] + real(yp + block_y) * vox_size[1];
		    real wy = vox_origin[1] + real(yp + block_y + y_step)
		      * vox_size[1];
		    sl_int xy_base = x_base + sl_int(block_y) * sl_int(nz);

		    for (int ax = 0; ax < a_step; ax++) {
		      int a = block_a + ax;
		      // rotate source and detector positions by current angle
		      real cphi = std::cos(angles[a]);
		      real sphi = std::sin(angles[a]);
	      
		      // from bproject
		      int hmin = 0;
		      int hmax = nh_pixels - 1;
		      if (std::abs(cphi) < epsilon) {
			if (sphi < 0.0) {
			  hmin = int(std::floor(vx * ihp_step - h_pix0));
			  if (hmin < 0)
			    hmin = 0;
			  hmax = int(std::floor((wx - epsilon) * ihp_step
						- h_pix0));
			  if (hmax >= nh_pixels)
			    hmax = nh_pixels - 1;
			} else {
			  hmin = int(std::ceil((- wx + epsilon) * ihp_step
					       - h_pix0));
			  if (hmin < 0)
			    hmin = 0;
			  hmax = int(std::ceil((- vx) * ihp_step - h_pix0));
			  if (hmax >= nh_pixels)
			    hmax = nh_pixels - 1;
			}
		      } else if (std::abs(sphi) < epsilon) {
			if (cphi < 0.0) {
			  hmin = int(std::ceil((- wy + epsilon) * ihp_step
					       - h_pix0));
			  if (hmin < 0)
			    hmin = 0;
			  hmax = int(std::ceil((- vy) * ihp_step - h_pix0));
			  if (hmax >= nh_pixels)
			    hmax = nh_pixels - 1;
			} else {
			  hmin = int(std::floor(vy * ihp_step - h_pix0));
			  if (hmin < 0)
			    hmin = 0;
			  hmax = int(std::floor((wy - epsilon) * ihp_step
						- h_pix0));
			  if (hmax >= nh_pixels)
			    hmax = nh_pixels - 1;
			}
		      } else {
			real ymin;
			real ymax;
			// Todo - pass in arrays of vy * cphi etc?
			if (cphi > 0.0) {
			  if (sphi > 0.0) {
			    ymin = vy * cphi - wx * sphi; //y10;
			    ymax = wy * cphi - vx * sphi; //y01;
			  } else {
			    ymin = vy * cphi - vx * sphi; //y00;
			    ymax = wy * cphi - wx * sphi; //y11;
			  }
			} else {
			  if (sphi > 0.0) {
			    ymin = wy * cphi - wx * sphi; //y11;
			    ymax = vy * cphi - vx * sphi; //y00;
			  } else {
			    ymin = wy * cphi - vx * sphi; //y01;
			    ymax = vy * cphi - wx * sphi; //y10;
			  }
			}
			if (ymin < h_pixels[0])
			  hmin = 0;
			else
			  hmin = int(std::floor(ymin * ihp_step - h_pix0));
			//if (ymax >= h_pixels[nh]) Todo
			hmax = std::min(int(std::floor(ymax * ihp_step
						       - h_pix0)),
					nh_pixels-1);
		      }
		      // end bproject
		
		      for (int h = hmin; h <= hmax; h++) {
			real p2_x = cphi * detector_x - sphi * h_pixels[h];
			real p2_y = sphi * detector_x + cphi * h_pixels[h];
			//real p1_x = p2_x - real(3.0) * cphi * detector_x;
			//real p1_y = p2_y - real(3.0) * sphi * detector_x;
			fproject_xy(p2_x, p2_y, pixels, voxels, vx, vy,
				    vox_size[0], vox_size[1], x_step, y_step,
				    nz, a, h, nv_pixels, d_conv, cphi, sphi,
				    ij_base, nyz, l_xy, ij_arr, ij_offsets,
				    block_yz, xy_base, nwork);
			if (nwork > 0) {
			  if (nwork > xy_size)
			    report_error("forward project overflow");
			  else {
			    machine::copy_to_device(ij_offsets.data(),
						    xy_offsets,
						    nwork * sizeof(int),
						    thread_id);
			    machine::copy_to_device(l_xy.data(), xy_buff,
						    nwork * sizeof(recon_type),
						    thread_id);
			    int pix_offset = (a * nh_pixels + h) * nv_pixels;
			    machine::run_parallel_xy(kernel_name, pix_buf,
						     pix_offset, vox_buf,
						     xy_buff, xy_offsets,
						     nwork, nv_pixels, nz,
						     nz, thread_id);
			    machine::accelerator_barrier(thread_id);
			  }
			}
		      }
		    }
		  }
		}
		machine::copy_from_device(pix_buf, &pixels[block_a][0][0],
					  a_step * sl_int(nh_pixels)
					  * sl_int(nv_pixels)
					  * sizeof(pixel_type), thread_id);
	      }
	    }
	  }
	}
	machine::device_free(xy_buff, thread_id);
	machine::device_free(xy_offsets, thread_id);
	machine::device_free(vox_buf, thread_id);
	machine::device_free(pix_buf, thread_id);
      }
    }
  }
}

void CCPi::parallel_beam::b2D_accel(const real_1d &h_pixels,
				    const real_1d &v_pixels,
				    const real_1d &angles, pixel_data &pixels,
				    voxel_data &voxels, const int n_angles,
				    const int nh_pixels, const int nv_pixels,
				    const real vox_origin[3],
				    const real vox_size[3],
				    const int nx, const int ny, const int nz)
{
  int_1d mapping(nv_pixels);
  int map_type = 0;
  gen_mapping(mapping, map_type, v_pixels, vox_origin[2], vox_size[2],
	      nv_pixels);
  if (map_type != 1 and map_type != 2)
    b2D_cpu(h_pixels, v_pixels, angles, pixels, voxels, n_angles, nh_pixels,
	    nv_pixels, vox_origin, vox_size, nx, ny, nz);
  else {
    // run a thread for each device
    int nthreads = machine::number_of_accelerators();
#pragma omp parallel shared(h_pixels, v_pixels, pixels, voxels, angles, vox_size, vox_origin, mapping) firstprivate(n_angles, nh_pixels, nv_pixels, nx, ny, nz, map_type, nthreads) num_threads(nthreads)
    {
      int thread_id = omp_get_thread_num();
      char kernel_name[32];
      if (map_type == 2)
	strcpy(kernel_name, "parallel_ah_z2");
      else
	strcpy(kernel_name, "parallel_ah_z1");
      if (!machine::check_accelerator_kernel(kernel_name, thread_id))
	report_error("Back kernel problem");
      else {
	real_1d c_angle(n_angles);
	real_1d s_angle(n_angles);
	for (int a = 0; a < n_angles; a++) {
	  real cos_phi = std::cos(angles[a]);
	  real sin_phi = std::sin(angles[a]);
	  c_angle[a] = cos_phi;
	  s_angle[a] = sin_phi;
	}

	recon_1d vox_z(nz + 1);
	for (int i = 0; i <= nz; i++)
	  vox_z[i] = vox_origin[2] + real(i) * vox_size[2];
	real_1d yvals(ny + 1);
	for (int j = 0; j <= ny; j++)
	  yvals[j] = vox_origin[1] + real(j) * vox_size[1];

	// rotation of voxel for offsets is same for all voxels
	real_1d y_offset(n_angles);
	real_1d i_offset(n_angles);
	real_1d length(n_angles);
	for (int a = 0; a < n_angles; a++) {
	  real cphi = c_angle[a];
	  real sphi = s_angle[a];
	  // use values for 0,0 voxel x_0 = b_x, x_1 = b_x + d_x etc
	  real y00 = vox_origin[1] * cphi - vox_origin[0] * sphi;
	  real y01 = (vox_origin[1] + vox_size[1]) * cphi - vox_origin[0] * sphi;
	  real y10 = vox_origin[1] * cphi - (vox_origin[0] + vox_size[0]) * sphi;
	  real y11 = (vox_origin[1] + vox_size[1]) * cphi
	    - (vox_origin[0] + vox_size[0]) * sphi;
	  real ymin;
	  real ybot;
	  if (cphi > 0.0) {
	    if (sphi > 0.0) {
	      ymin = y10;
	      if (y11 < y00) {
		ybot = y11;
	      } else {
		ybot = y00;
	      }
	    } else {
	      ymin = y00;
	      if (y01 < y10) {
		ybot = y01;
	      } else {
		ybot = y10;
	      }
	    }
	  } else {
	    if (sphi > 0.0) {
	      ymin = y11;
	      if (y01 < y10) {
		ybot = y01;
	      } else {
		ybot = y10;
	      }
	    } else {
	      ymin = y01;
	      if (y11 < y00) {
		ybot = y11;
	      } else {
		ybot = y00;
	      }
	    }
	  }
	  y_offset[a] = ybot - ymin;
	  i_offset[a] = 1.0 / y_offset[a];
	  if (std::abs(cphi) > std::abs(sphi))
	    length[a] = vox_size[0] / std::abs(cphi);
	  else
	    length[a] = vox_size[1] / std::abs(sphi);
	}
	const real ihp_step = 1.0 / (h_pixels[1] - h_pixels[0]);
	const real h_pix0 = h_pixels[0] / (h_pixels[1] - h_pixels[0]);

	// calc space on accelerator for copy
	sl_int accel_size =
	  machine::largest_alloc(thread_id) / sizeof(pixel_type);
	sl_int accel_proj = accel_size / (nh_pixels * nv_pixels);
	if (accel_proj > n_angles)
	  accel_proj = n_angles;

	const int x_block = 32;
	const int y_block = 32;
	const int a_block = accel_proj;

	// from bproject_ah
	const int pix_per_vox = nv_pixels / (nz - 1);
	// How big should the array be - Todo - use mapping for pix_per_vox?
	int nwork = 0;
	int xy_size = 2 * pix_per_vox * (a_block + 10);
	if (xy_size % 32 != 0)
	  xy_size += 32 - (xy_size % 32);
	pixel_ptr_1d ah_arr(xy_size);
	int_1d ah_offsets(xy_size);
	recon_1d l_xy(xy_size);
	// Need 3D as can't overwrite 1D ones while doing async copy.
	boost::multi_array<int, 3> ah3_offsets(boost::extents[x_block][y_block][xy_size]);
	boost::multi_array<recon_type, 3> l3_xy(boost::extents[x_block][y_block][xy_size]);
	boost::multi_array<int, 2> work_sizes(boost::extents[x_block][y_block]);

	dev_ptr pix_buf = machine::device_allocate(accel_proj
						   * sl_int(nh_pixels)
						   * sl_int(nv_pixels)
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
	dev_ptr xy_buff = machine::device_allocate(sl_int(x_block)
						   * sl_int(y_block) * xy_size
						   * sizeof(recon_type),
						   true, thread_id);
	dev_ptr xy_work = machine::device_allocate(sl_int(x_block)
						   * sl_int(y_block)
						   * sizeof(int),
						   true, thread_id);

	long counter = 0;
	std::vector<event_t> vox_x_ev(x_block);
	for (int ap = 0; ap < n_angles; ap += accel_proj) {
	  int p_step = accel_proj;
	  if (ap + p_step > n_angles)
	    p_step = n_angles - ap;
	  // Todo - should be contiguous as block of ap
	  event_t pixel_ev;
	  machine::copy_to_device(&pixels[ap][0][0], pix_buf,
				  p_step * sl_int(nh_pixels) * sl_int(nv_pixels)
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
		// Its not contiguous, since not a full y range
		for (int ix = 0; ix < x_step; ix++)
		  machine::copy_to_device(&voxels[block_x + ix][block_y][0],
					  vox_buf, sl_int(ix * y_step)
					  * sl_int(nz) * sizeof(voxel_type),
					  sl_int(y_step)
					  * sl_int(nz) * sizeof(voxel_type),
					  thread_id, &vox_x_ev[ix]);
		for (int block_a = 0; block_a < p_step; block_a += a_block) {
		  int a_step = a_block;
		  if (block_a + a_step > p_step)
		    a_step = p_step - block_a;
		  // we don't block h since we have a min/max from the
		  // x/y blocks
		  int offset = 0;
		  std::vector<std::vector<event_t> *> ev(x_block);
		  for (int ix = 0; ix < x_step; ix++) {
		    ev[ix] = new std::vector<event_t>(5);
		    int i = block_x + ix;
		    const real x_0 = vox_origin[0] + real(i) * vox_size[0];
		    const real x_n = vox_origin[0] + real(i + 1) * vox_size[0];
		    for (int jx = 0; jx < y_step; jx++) {
		      int j = block_y + jx;
		      bproject_ah(pixels, voxels, x_0, yvals[j],
				  x_n, yvals[j + 1], vox_size[0], vox_size[1],
				  nz, i, j, a_step, nh_pixels, nv_pixels,
				  h_pixels, c_angle, s_angle, y_offset,
				  i_offset, length, h_pix0, ihp_step,
				  ap + block_a, ah_arr, l_xy, ah_offsets,
				  block_a, nwork);
		      offset += xy_size;
		      if (nwork > 0) {
			if (nwork > xy_size) {
			  report_error("back project overflow");
			  nwork = 0;
			} else {
			  for (int d = 0; d < nwork; d++)
			    ah3_offsets[ix][jx][d] = ah_offsets[d];
			  for (int d = 0; d < nwork; d++)
			    l3_xy[ix][jx][d] = l_xy[d];
			}
		      }
		      work_sizes[ix][jx] = nwork;
		    }
		    // Todo - can we use a local event array or do we need
		    // to keep each one while the async run might occur?
		    // e.g. does opencl have its own copy
		    machine::copy_to_device(&ah3_offsets[ix][0][0], xy_offsets,
					    ix * y_step * xy_size * sizeof(int),
					    y_step * xy_size * sizeof(int),
					    thread_id, &((*(ev[ix]))[0]));
		    machine::copy_to_device(&l3_xy[ix][0][0], xy_buff,
					    ix * y_step * xy_size
					    * sizeof(recon_type),
					    y_step * xy_size
					    * sizeof(recon_type),
					    thread_id, &((*(ev[ix]))[1]));
		    machine::copy_to_device(&work_sizes[ix][0], xy_work,
					    ix * y_step * sizeof(int),
					    y_step * sizeof(int), thread_id,
					    &((*(ev[ix]))[2]));
		    (*(ev[ix]))[3] = pixel_ev;
		    (*(ev[ix]))[4] = vox_x_ev[ix];
		    int vox_offset = (ix * y_step) * nz;
		    machine::run_parallel_ah(kernel_name, pix_buf,
					     vox_buf, vox_offset, xy_buff,
					     xy_offsets, xy_work, nv_pixels,
					     nz, xy_size, ix, nz, y_step,
					     thread_id, ev[ix]);
		  }
		  machine::accelerator_barrier(thread_id);
		  for (int ix = 0; ix < x_step; ix++)
		    delete ev[ix];
		}
		// no events copy is blocking
		for (int ix = 0; ix < x_step; ix++)
		  machine::copy_from_device(vox_buf,
					    &voxels[block_x + ix][block_y][0],
					    sl_int(ix * y_step)
					    * sl_int(nz) * sizeof(voxel_type),
					    sl_int(y_step)
					    * sl_int(nz) * sizeof(voxel_type),
					    thread_id);
		//machine::accelerator_barrier(thread_id);
	      }
	      counter++;
	    }
	  }
	}
	machine::device_free(xy_buff, thread_id);
	machine::device_free(xy_offsets, thread_id);
	machine::device_free(vox_buf, thread_id);
	machine::device_free(pix_buf, thread_id);
      }
    }
  }
}

#else

void CCPi::parallel_beam::f2D_accel(const real_1d &h_pixels,
				    const real_1d &v_pixels,
				    const real_1d &angles, const int n_angles,
				    const int nh_pixels, const int nv_pixels,
				    const real vox_origin[3],
				    const real vox_size[3],
				    const int nx, const int ny, const int nz,
				    pixel_data &pixels, voxel_data &voxels)
{
  report_error("BUG: attempt to run on accelerator");
}

void CCPi::parallel_beam::b2D_accel(const real_1d &h_pixels,
				    const real_1d &v_pixels,
				    const real_1d &angles, pixel_data &pixels,
				    voxel_data &voxels, const int n_angles,
				    const int nh_pixels, const int nv_pixels,
				    const real vox_origin[3],
				    const real vox_size[3],
				    const int nx, const int ny, const int nz)
{
  report_error("BUG: attempt to run on accelerator");
}

#endif // USE_OPENCL
