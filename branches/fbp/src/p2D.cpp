
#include <map>
#include <vector>
#include <cmath>
#include <float.h>
#ifdef MATLAB_MEX_FILE
#  include "mex_types.hpp"
#else
#  include "base_types.hpp"
#endif // mex
#include "filters.hpp"
#include "instruments.hpp"
#include "ui_calls.hpp"
#ifdef TEST2D
#  include <iostream>
#endif // TEST2D
#if defined(__AVX__)
#  include <immintrin.h>
#endif // AVX etc

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
  const recon_type *lptr = assume_aligned(&(l_xy[0]), recon_type);
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
      // high to low - abcdabcd->aabbccdd
      const __m256i offsets = _mm256_set_epi32(3, 3, 2, 2, 1, 1, 0, 0);
      for (int m = 0; m < n; m++) {
	const voxel_type *const vox = voxels[m];
	// __mm256 al = _mm256_broadcast_ss(lptr[m]); ?
	const recon_type alpha = lptr[m];
	__m256 al = _mm256_set1_ps(alpha);
	int h = 0;
	for (int l = 0; l < nv; l += 8) {
	  __m256 vx = _mm256_broadcast_ps((__m128 *)&vox[h]); //abcdabcd
	  __m256 v = _mm256_permutevar_ps(vx, offsets);
	  _mm256_store_ps(&pix[l], _mm256_fmadd_ps(al, v,
						   _mm256_load_ps(&pix[l])));
	  h += 4;
	}
      }
    }
#elif defined(__AVX__) && PIXEL_SIZE == 4
    {
      const __m256i offsets = _mm256_set_epi32(3, 3, 2, 2, 1, 1, 0, 0);
      for (int m = 0; m < n; m++) {
	const voxel_type *const vox = voxels[m];
	// __mm256 al = _mm256_broadcast_ss(lptr[m]); ?
	const recon_type alpha = lptr[m];
	__m256 al = _mm256_set1_ps(alpha);
	int h = 0;
	// Todo - can load 8 way vox as 2x means 16 way pix unroll
	for (int l = 0; l < nv; l += 8) {
	  __m256 vx = _mm256_broadcast_ps((__m128 *)&vox[h]); //abcdabcd
	  __m256 v = _mm256_permutevar_ps(vx, offsets);
	  _mm256_store_ps(&pix[l], _mm256_add_ps(_mm256_mul_ps(al, v),
						 _mm256_load_ps(&pix[l])));
	  h += 4;
	}
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
				      const sl_int nyz, const int_1d &mapping,
				      const int map_type)
{
  // Todo? In parallel beam some of this should be common in a since
  // all h within a have same angle to voxels.
  int max_n = std::max(nx, ny);
  recon_1d l_xy(2 * max_n);
  voxel_ptr_1d ij_arr(2 * max_n + 1);
  int count = 0;
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
	  for (int j = ny - 1; j >= 0; j--) {
	    l_xy[count] = d_y;
	    ij_arr[count] = ij_offset;
	    ij_offset -= nz;
	    count++;
	  }
	} else {
	  voxel_ptr ij_offset = voxels.data() + ij_base + sl_int(i) * nyz;
	  for (int j = 0; j < ny; j++) {
	    l_xy[count] = d_y;
	    ij_arr[count] = ij_offset;
	    ij_offset += nz;
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
	for (int i = nx - 1; i >= 0; i--) {
	  l_xy[count] = d_x;
	  ij_arr[count] = ij_offset;
	  ij_offset -= nyz;
	  count++;
	}
      } else {
	voxel_ptr ij_offset = voxels.data() + ij_base + sl_int(j) * sl_int(nz);
	for (int i = 0; i < nx; i++) {
	  l_xy[count] = d_x;
	  ij_arr[count] = ij_offset;
	  ij_offset += nyz;
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
	  voxel_ptr xy_offset = voxels.data() + ij_base + sl_int(x) * nyz
	    + sl_int(y) * sl_int(nz);
	  while (x < nx and y < ny) {
	    ij_arr[count] = xy_offset;
	    if (alpha_x[x + 1] < alpha_y[y + 1] - epsilon) {
	      l_xy[count] = (alpha_x[x + 1] - alpha_p);
	      alpha_p = alpha_x[x + 1];
	      xy_offset += nyz;
	      x++;
	    } else if (alpha_y[y + 1] < alpha_x[x + 1] - epsilon) {
	      l_xy[count] = (alpha_y[y + 1] - alpha_p);
	      alpha_p = alpha_y[y + 1];
	      xy_offset += nz;
	      y++;
	    } else {
	      real mx = std::max(alpha_x[x + 1], alpha_y[y + 1]);
	      l_xy[count] = (mx - alpha_p);
	      x++;
	      y++;
	      xy_offset += nyz + nz;
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
	  voxel_ptr xy_offset = voxels.data() + ij_base + sl_int(x) * nyz
	    + sl_int(y) * sl_int(nz);
	  while (x < nx and y >= 0) {
	    ij_arr[count] = xy_offset;
	    if (alpha_x[x + 1] < alpha_y[y] - epsilon) {
	      l_xy[count] = (alpha_x[x + 1] - alpha_p);
	      alpha_p = alpha_x[x + 1];
	      xy_offset += nyz;
	      x++;
	    } else if (alpha_y[y] < alpha_x[x + 1] - epsilon) {
	      l_xy[count] = (alpha_y[y] - alpha_p);
	      alpha_p = alpha_y[y];
	      xy_offset -= nz;
	      y--;
	    } else {
	      real mx = std::max(alpha_x[x + 1], alpha_y[y]);
	      l_xy[count] = (mx - alpha_p);
	      x++;
	      y--;
	      xy_offset += nyz - nz;
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
	  voxel_ptr xy_offset = voxels.data() + ij_base + sl_int(x) * nyz
	    + sl_int(y) * sl_int(nz);
	  while (x >= 0 and y < ny) {
	    ij_arr[count] = xy_offset;
	    if (alpha_x[x] < alpha_y[y + 1] - epsilon) {
	      l_xy[count] = (alpha_x[x] - alpha_p);
	      alpha_p = alpha_x[x];
	      xy_offset -= nyz;
	      x--;
	    } else if (alpha_y[y + 1] < alpha_x[x] - epsilon) {
	      l_xy[count] = (alpha_y[y + 1] - alpha_p);
	      alpha_p = alpha_y[y + 1];
	      xy_offset += nz;
	      y++;
	    } else {
	      real mx = std::max(alpha_y[y + 1], alpha_x[x]);
	      l_xy[count] = (mx - alpha_p);
	      x--;
	      y++;
	      xy_offset += -nyz + nz;
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
	  voxel_ptr xy_offset = voxels.data() + ij_base + sl_int(x) * nyz
	    + sl_int(y) * sl_int(nz);
	  while (x >= 0 and y >= 0) {
	    ij_arr[count] = xy_offset;
	    if (alpha_x[x] < alpha_y[y] - epsilon) {
	      l_xy[count] = (alpha_x[x] - alpha_p);
	      alpha_p = alpha_x[x];
	      xy_offset -= nyz;
	      x--;
	    } else if (alpha_y[y] < alpha_x[x] - epsilon) {
	      l_xy[count] = (alpha_y[y] - alpha_p);
	      alpha_p = alpha_y[y];
	      xy_offset -= nz;
	      y--;
	    } else {
	      real mx = std::max(alpha_x[x], alpha_y[y]);
	      l_xy[count] = (mx - alpha_p);
	      x--;
	      y--;
	      xy_offset -= (nyz + nz);
	      alpha_p = mx;
	    }
	    count++;
	  }
	}
      }
    }
  }
  if (count > 2 * max_n + 1)
    report_error("forward project overflow");
  if (count > 0) {
    calc_xy_z(&(pixels[a][h][0]), ij_arr, l_xy, count, nv, nz,
	      mapping, map_type);
  }    
}

void CCPi::parallel_beam::f2D(const real_1d &h_pixels, const real_1d &v_pixels,
			      const real_1d &angles, const int n_angles,
			      const int nh_pixels, const int nv_pixels,
			      const real vox_origin[3], const real vox_size[3],
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
			d_conv, cphi, sphi, ij_base, nyz, mapping, map_type);
	  }
	}
      }
    }
  }
}

void CCPi::parallel_beam::calc_ah_z(const pixel_ptr_1d &pixels,
				    voxel_type *const voxels,
				    const recon_1d &l_xy, const int n,
				    const int nv, const int nz,
				    const int_1d &mapping, const int map_type)
{
  voxel_type *const vox = assume_aligned(voxels, voxel_type);
  const recon_type *lptr = assume_aligned(&(l_xy[0]), recon_type);
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
	__m256 al = _mm256_set1_ps(alpha);
	int v = 0;
	for (int l = 0; l < nz; l += 8) {
	  __m256 pl = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_load_ps(&pix[v + 0])),
					   _mm_load_ps(&pix[v + 8]), 0x1);
	  __m256 ph = _mm256_insertf128_ps(_mm256_castps128_ps256(_mm_load_ps(&pix[v + 4])),
					   _mm_load_ps(&pix[v + 12]), 0x1);
	  _mm256_store_ps(&vox[l], _mm256_fmadd_ps(al, _mm256_hadd_ps(pl, ph),
						   _mm256_load_ps(&vox[l])));
	  v += 16;
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
	__m256 pl = _mm256_hadd_ps(pa, pb);
	__m256 ph = _mm256_hadd_ps(pc, pd);
	_mm256_store_ps(&vox[l], _mm256_fmadd_ps(al, _mm256_hadd_ps(pl, ph),
						 _mm256_load_ps(&vox[l])));
	v += 32;
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
				      const int_1d &mapping, const int map_type)
{
  // Rather than using the centre just calculate for all 4 corners,
  // generate h values and loop from smallest to largest.
  const int pix_per_vox = n_v / (nz - 1);
  // How big should the array be - Todo - use mapping for pix_per_vox?
  int count = 0;
  pixel_ptr_1d ah_arr(2 * pix_per_vox * (n_angles + 10));
  recon_1d l_xy(2 * pix_per_vox * (n_angles + 10));
  // corners (x0,y0), (x0,yn), (xn,y0), (xn,yn)
  // Todo - in parallel we can probably make a better guess at which 2 corners
  // we need for the upper and lower limits.
  //real pixel_step = h_pixels[1] - h_pixels[0];
  sl_int nah = sl_int(n_h) * sl_int(n_v);
  pixel_ptr ah_offset = pixels.data() + sl_int(a_off) * nah;
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
	for (int h = hmin; h <= hmax; h++) {
	  if (h_pixels[h] >= x_0 and h_pixels[h] < x_n) {
	    l_xy[count] = d_y;
	    ah_arr[count] = h_offset;
	    count++;
	  }
	  h_offset += n_v;
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
	for (int h = hmin; h <= hmax; h++) {
	  if (-h_pixels[h] >= x_0 and -h_pixels[h] < x_n) {
	    l_xy[count] = d_y;
	    ah_arr[count] = h_offset;
	    count++;
	  }
	  h_offset += n_v;
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
	for (int h = hmin; h <= hmax; h++) {
	  if (- h_pixels[h] >= y_0 and - h_pixels[h] < y_n) {
	    l_xy[count] = d_x;
	    ah_arr[count] = h_offset;
	    count++;
	  }
	  h_offset += n_v;
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
	for (int h = hmin; h <= hmax; h++) {
	  if (h_pixels[h] >= y_0 and h_pixels[h] < y_n) {
	    l_xy[count] = d_x;
	    ah_arr[count] = h_offset;
	    count++;
	  }
	  h_offset += n_v;
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
      for (int h = hmin; h <= hmax; h++) {
	if (h_pixels[h] > ymin + epsilon) {
	  if (h_pixels[h] < ybot) {
	    l_xy[count] = l * (h_pixels[h] - ymin) * i_offset[a];
	    ah_arr[count] = h_offset;
	    count++;
	  } else if (h_pixels[h] <= ytop) {
	    l_xy[count] = l;
	    ah_arr[count] = h_offset;
	    count++;
	  } else if (h_pixels[h] < ymax - epsilon) {
	    l_xy[count] = l * (ymax - h_pixels[h]) * i_offset[a];
	    ah_arr[count] = h_offset;
	    count++;
	  } else
	    break;
	}
	h_offset += n_v;
      }
    }
    ah_offset += nah;
  }
  if (count > 2 * pix_per_vox * (n_angles + 10))
    report_error("back project overflow");
  if (count > 0) {
    calc_ah_z(ah_arr, &(voxels[i][j][0]), l_xy, count, n_v, nz,
	      mapping, map_type);
  }
}

void CCPi::parallel_beam::b2D(const real_1d &h_pixels, const real_1d &v_pixels,
			      const real_1d &angles, pixel_data &pixels,
			      voxel_data &voxels, const int n_angles,
			      const int nh_pixels, const int nv_pixels,
			      const real vox_origin[3], const real vox_size[3],
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
			h_pix0, ihp_step, block_a, mapping, map_type);
	  }
	}
      }
    }
  }      
}
