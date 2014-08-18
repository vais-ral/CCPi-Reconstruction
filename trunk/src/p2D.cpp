
#include <map>
#include <vector>
#include <cmath>
#include <float.h>
#ifdef MATLAB_MEX_FILE
#  include "mex_types.hpp"
#else
#  include "base_types.hpp"
#endif // mex
#include "instruments.hpp"
#include "ui_calls.hpp"
#ifdef TEST2D
#  include <iostream>
#endif // TEST2D

static const recon_type epsilon = FLT_EPSILON;

extern bool test_2D(const real start[], const real end[],
		    const real b_x, const real b_y,
		    const real d_x, const real d_y,
		    const int im_size_x, const int im_size_y,
		    recon_type &length1, recon_type &length2);
extern bool test_3D(const real start[], const real end[],
		    const real b_x, const real b_y, const real b_z,
		    const real d_x, const real d_y, const real d_z,
		    const int im_size_x, const int im_size_y,
		    const int im_size_z, recon_type &length);

void CCPi::parallel_beam::gen_mapping(std::vector<int> &mapping, int &map_type,
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

void CCPi::parallel_beam::calc_xy_z(pixel_data &pixels, voxel_data &voxels,
				    const recon_1d &l_xy,
				    const std::vector<long> &ij,
				    const int n, const int a, const int h,
				    const int nv, const int nz,
				    const std::vector<int> &mapping,
				    const int map_type)
{
  pixel_type *const pix = &(pixels[a][h][0]);
  switch (map_type) {
  case 1:
    for (int m = 0; m < n; m++) {
      const voxel_type *const vox = &(voxels.data()[ij[m]]);
      const recon_type alpha = l_xy[m];
      for (int v = 0; v < nv; v++)
	pix[v] += vox[v] * alpha;
    }
    break;
  case 2:
    {
      int v2 = nv / 2;
      for (int m = 0; m < n; m++) {
	const voxel_type *const vox = &(voxels.data()[ij[m]]);
	const recon_type alpha = l_xy[m];
	int v = 0;
	for (int l = 0; l < v2; l++) {
	  pix[v + 0] += vox[l] * alpha;
	  pix[v + 1] += vox[l] * alpha;
	  v += 2;
	}
      }
    }
    break;
  case 4:
    {
      int v4 = nv / 4;
      for (int m = 0; m < n; m++) {
	const voxel_type *const vox = &(voxels.data()[ij[m]]);
	const recon_type alpha = l_xy[m];
	int v = 0;
	for (int l = 0; l < v4; l++) {
	  pix[v + 0] += vox[l] * alpha;
	  pix[v + 1] += vox[l] * alpha;
	  pix[v + 2] += vox[l] * alpha;
	  pix[v + 3] += vox[l] * alpha;
	  v += 4;
	}
      }
    }
    break;
  default:
    for (int m = 0; m < n; m++) {
      const voxel_type *const vox = &(voxels.data()[ij[m]]);
      const recon_type alpha = l_xy[m];
      for (int v = 0; v < nv; v++)
	pix[v] += vox[mapping[v]] * alpha;
    }
    break;
  }
}

void CCPi::parallel_beam::fproject_xy(const real p1_x, const real p1_y,
				      const real p2_x, const real p2_y,
				      pixel_data &pixels, voxel_data &voxels,
				      const real b_x, const real b_y,
				      const real b_z, const real d_x,
				      const real d_y, const real d_z,
				      const int nx, const int ny, const int nz,
				      const int a, const int h,
				      const int nv, const recon_type d_conv,
				      const real_1d &v_pixels,
				      const real cphi, const real sphi,
				      const std::vector<int> &mapping,
				      const int map_type)
{
  // Todo? In parallel beam some of this should be common in a since
  // all h within a have same angle to voxels.
  int max_n = std::max(nx, ny);
  recon_1d l_xy(2 * max_n);
  std::vector<long> ij_arr(2 * max_n + 1);
  int count = 0;
  long nyz = long(ny) * long(nz);
  if (std::abs(cphi) < epsilon) {
    if (std::abs(sphi) < epsilon) {
      // Its not a line - shouldn't happen
      //return;
    } else {
      // line parallel to y
      int i = int(std::floor((p2_x - b_x) / d_x)); // == q2_x
      if (i >= 0 and i < nx) {
	if (sphi < 0.0) {
	  long ij_offset = long(i) * nyz + long(ny - 1) * long(nz);
	  for (int j = ny - 1; j >= 0; j--) {
	    l_xy[count] = d_y;
	    ij_arr[count] = ij_offset;
	    ij_offset -= nz;
	    count++;
	  }
	} else {
	  long ij_offset = long(i) * nyz;
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
	long ij_offset = long(nx - 1) * nyz + long(j) * long(nz);
	for (int i = nx - 1; i >= 0; i--) {
	  l_xy[count] = d_x;
	  ij_arr[count] = ij_offset;
	  ij_offset -= nyz;
	  count++;
	}
      } else {
	long ij_offset = long(j) * long(nz);
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
	  long xy_offset = long(x) * nyz + long(y) * long(nz);
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
	  long xy_offset = long(x) * nyz + long(y) * long(nz);
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
	  long xy_offset = long(x) * nyz + long(y) * long(nz);
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
	  long xy_offset = long(x) * nyz + long(y) * long(nz);
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
#ifdef TEST2D
  {
    std::vector<bool> cmp(count);
    for (int i = 0; i < count; i++)
      cmp[i] = false;
    real start[3];
    real end[3];
    start[0] = p1_x;
    start[1] = p1_y;
    end[0] = p2_x;
    end[1] = p2_y;
    for (int i = 0; i < nx; i++) {
      real x_0 = b_x + real(i) * d_x;
      for (int j = 0; j < ny; j++) {
	real y_0 = b_y + real(j) * d_y;
	recon_type ln1, ln2;
	if (test_2D(start, end, x_0, y_0, d_x, d_y, 1, 1, ln1, ln2)) {
	  recon_type ln = ln1 - ln2;
	  int k;
	  for (k = 0; k < count; k++) {
	    if (ij_arr[k] == long(i) * nyz + long(j) * long(nz)) {
	      cmp[k] = true;
	      real diff = l_xy[k] / d_conv;
	      if (ln < diff - epsilon or ln > diff + epsilon)
		std::cerr << "Bug " << i << ' ' << j << ' ' << a << ' ' << h 
			  << ' ' << ln << ' ' << diff << '\n';
	      break;
	    }
	  }
	  if (k == count and ln > epsilon)
	    std::cerr << "Missed " << i << ' ' << j << ' ' << a << ' ' << h
		      << ' ' << k << ' ' << count << ' ' << ln << '\n';
	}
      }
    }
    for (int k = 0; k < count; k++) {
      if (!cmp[k])
	std::cerr << "Not found " << a << ' ' << h << ' ' << ij_arr[k]
		  << ' ' << ij_arr[k] << ' ' << k << ' ' << count
		  << ' ' << l_xy[k] / d_conv << '\n';
    }
  }
#endif // TEST2D
  if (count > 0) {
    calc_xy_z(pixels, voxels, l_xy, ij_arr, count, a, h,
	      nv, nz, mapping, map_type);
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

  std::vector<int> mapping(nv_pixels);
  int map_type = 0;
  gen_mapping(mapping, map_type, v_pixels, vox_origin[2], vox_size[2],
	      nv_pixels);
  
#pragma omp parallel for shared(h_pixels, v_pixels, pixels, voxels, angles, d_conv, vox_size, vox_origin, mapping) firstprivate(n_angles, nh_pixels, nv_pixels, nx, ny, nz, detector_x, map_type) schedule(dynamic)
  for (int a = 0; a < n_angles; a++) {
    for (int h = 0; h < nh_pixels; h++) {
      // rotate source and detector positions by current angle
      real cosa = std::cos(angles[a]);
      real sina = std::sin(angles[a]);

      real p2_x = cosa * detector_x - sina * h_pixels[h];
      real p2_y = sina * detector_x + cosa * h_pixels[h];
      real p1_x = p2_x - real(3.0) * cosa * detector_x;
      real p1_y = p2_y - real(3.0) * sina * detector_x;
      fproject_xy(p1_x, p1_y, p2_x, p2_y, pixels, voxels, vox_origin[0],
		  vox_origin[1], vox_origin[2], vox_size[0], vox_size[1],
		  vox_size[2], nx, ny, nz, a, h, nv_pixels, d_conv, v_pixels,
		  cosa, sina, mapping, map_type);
    }
  }
}

void CCPi::parallel_beam::calc_ah_z(pixel_data &pixels, voxel_data &voxels,
				    const recon_1d &l_xy,
				    const std::vector<long> &ah,
				    const int n, const int i, const int j,
				    const int nv, const int nz,
				    const std::vector<int> &mapping,
				    const int map_type)
{
  voxel_type *const vox = &(voxels[i][j][0]);
  switch (map_type) {
  case 1:
    for (int m = 0; m < n; m++) {
      const pixel_type *const pix = &(pixels.data()[ah[m]]);
      const recon_type alpha = l_xy[m];
      for (int v = 0; v < nv; v++)
	vox[v] += pix[v] * alpha;
    }
    break;
  case 2:
    {
      int v2 = nv / 2;
      for (int m = 0; m < n; m++) {
	const pixel_type *const pix = &(pixels.data()[ah[m]]);
	const recon_type alpha = l_xy[m];
	int v = 0;
	for (int l = 0; l < v2; l++) {
	  vox[l] += pix[v + 0] * alpha;
	  vox[l] += pix[v + 1] * alpha;
	  v += 2;
	}
      }
    }
    break;
  case 4:
    {
      int v4 = nv / 4;
      for (int m = 0; m < n; m++) {
	const pixel_type *const pix = &(pixels.data()[ah[m]]);
	const recon_type alpha = l_xy[m];
	int v = 0;
	for (int l = 0; l < v4; l++) {
	  vox[l] += pix[v + 0] * alpha;
	  vox[l] += pix[v + 1] * alpha;
	  vox[l] += pix[v + 2] * alpha;
	  vox[l] += pix[v + 3] * alpha;
	  v += 4;
	}
      }
    }
    break;
  default:
    for (int m = 0; m < n; m++) {
      const pixel_type *const pix = &(pixels.data()[ah[m]]);
      const recon_type alpha = l_xy[m];
      for (int v = 0; v < nv; v++)
	vox[mapping[v]] += pix[v] * alpha;
    }
    break;
  }
}

void CCPi::parallel_beam::bproject_ah(const real source_x,
				      const real detector_x, pixel_data &pixels,
				      voxel_data &voxels, const real x_0,
				      const real y_0, const real x_n,
				      const real y_n, const real b_z,
				      const real d_x, const real d_y,
				      const real d_z, const int nx,
				      const int ny, const int nz, const int i,
				      const int j, const int n_angles,
				      const int n_h, const int n_v,
				      const real_1d &h_pixels,
				      const real_1d &v_pixels,
				      const real_1d &cangle,
				      const real_1d &sangle,
				      const real_1d &y_offset,
				      const real_1d &i_offset,
				      const real_1d &length,
				      const real h_pix0, const real ihp_step,
				      const recon_type d_conv,
				      const std::vector<int> &mapping,
				      const int map_type)
{
  // Rather than using the centre just calculate for all 4 corners,
  // generate h values and loop from smallest to largest.
  const int pix_per_vox = n_v / (nz - 1);
  // How big should the array be - Todo - use mapping for pix_per_vox?
  int count = 0;
  std::vector<long> ah_arr(2 * pix_per_vox * n_angles);
  recon_1d l_xy(2 * pix_per_vox * n_angles);
  // corners (x0,y0), (x0,yn), (xn,y0), (xn,yn)
  // Todo - in parallel we can probably make a better guess at which 2 corners
  // we need for the upper and lower limits.
  //real pixel_step = h_pixels[1] - h_pixels[0];
  long nah = long(n_h) * long(n_v);
  long ah_offset = 0;
  for (int a = 0; a < n_angles; a++) {
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
	long h_offset = long(hmin) * long(n_v); 
	for (int h = hmin; h <= hmax; h++) {
	  if (h_pixels[h] >= x_0 and h_pixels[h] < x_n) {
	    l_xy[count] = d_y;
	    ah_arr[count] = ah_offset + h_offset;
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
	long h_offset = long(hmin) * long(n_v); 
	for (int h = hmin; h <= hmax; h++) {
	  if (-h_pixels[h] >= x_0 and -h_pixels[h] < x_n) {
	    l_xy[count] = d_y;
	    ah_arr[count] = ah_offset + h_offset;
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
	long h_offset = long(hmin) * long(n_v); 
	for (int h = hmin; h <= hmax; h++) {
	  if (- h_pixels[h] >= y_0 and - h_pixels[h] < y_n) {
	    l_xy[count] = d_x;
	    ah_arr[count] = ah_offset + h_offset;
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
	long h_offset = long(hmin) * long(n_v); 
	for (int h = hmin; h <= hmax; h++) {
	  if (h_pixels[h] >= y_0 and h_pixels[h] < y_n) {
	    l_xy[count] = d_x;
	    ah_arr[count] = ah_offset + h_offset;
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
      long h_offset = long(hmin) * long(n_v);
      for (int h = hmin; h <= hmax; h++) {
	if (h_pixels[h] > ymin + epsilon) {
	  if (h_pixels[h] < ybot) {
	    l_xy[count] = l * (h_pixels[h] - ymin) * i_offset[a];
	    ah_arr[count] = ah_offset + h_offset;
	    count++;
	  } else if (h_pixels[h] <= ytop) {
	    l_xy[count] = l;
	    ah_arr[count] = ah_offset + h_offset;
	    count++;
	  } else if (h_pixels[h] < ymax - epsilon) {
	    l_xy[count] = l * (ymax - h_pixels[h]) * i_offset[a];
	    ah_arr[count] = ah_offset + h_offset;
	    count++;
	  } else
	    break;
	}
	h_offset += n_v;
      }
    }
    ah_offset += nah;
  }
  if (count > 2 * pix_per_vox * n_angles)
    report_error("back project overflow");
#ifdef TEST2D
  {
    std::vector<bool> cmp(count);
    for (int k = 0; k < count; k++)
      cmp[k] = false;
    real start[3];
    real end[3];
    for (int a = 0; a < n_angles; a++) {
      for (int h = 0; h < n_h; h++) {
	real cos_curr_angle = cangle[a];
	real sin_curr_angle = sangle[a];
	end[0] = cos_curr_angle * detector_x - sin_curr_angle * h_pixels[h];
	end[1] = sin_curr_angle * detector_x + cos_curr_angle * h_pixels[h];
	start[0] = end[0] - real(3.0) * cos_curr_angle * detector_x;
	start[1] = end[1] - real(3.0) * sin_curr_angle * detector_x;
	recon_type ln1, ln2;
	if (test_2D(start, end, x_0, y_0, d_x, d_y, 1, 1, ln1, ln2)) {
	  recon_type ln = ln1 - ln2;
	  int k;
	  for (k = 0; k < count; k++) {
	    if (ah_arr[k] == long(a) * nah + long(h) * long(n_v)) {
	      cmp[k] = true;
	      real diff = l_xy[k] / d_conv;
	      if (ln < diff - epsilon or ln > diff + epsilon)
		std::cerr << "Bug " << a << ' ' << h << ' ' << i << ' ' << j 
			  << ' ' << ln << ' ' << diff << '\n';
	      break;
	    }
	  }
	  if (k == count and ln > epsilon)
	    std::cerr << "Missed " << a << ' ' << h << ' ' << i << ' ' << j
		      << ' ' << k << ' ' << count << ' ' << ln << '\n';
	}
      }
    }
    for (int k = 0; k < count; k++) {
      if (!cmp[k])
	std::cerr << "Not found " << i << ' ' << j << ' ' << ah_arr[k]
		  << ' ' << ah_arr[k] << ' ' << k << ' ' << count
		  << ' ' << l_xy[k] / d_conv << '\n';
    }
  }
#endif // TEST2D
  if (count > 0) {
    calc_ah_z(pixels, voxels, l_xy, ah_arr, count,
	      i, j, n_v, nz, mapping, map_type);
  }
}

void CCPi::parallel_beam::b2D(const real_1d &h_pixels,
			      const real_1d &v_pixels,
			      const real_1d &angles,
			      pixel_data &pixels,
			      voxel_data &voxels,
			      const int n_angles,
			      const int nh_pixels,
			      const int nv_pixels,
			      const real vox_origin[3],
			      const real vox_size[3],
			      const int nx,
			      const int ny,
			      const int nz)
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
  real detector_x = real(2.0) * std::max(std::abs(vox_origin[0]),
					 std::max(std::abs(vox_origin[1]),
						  std::abs(vox_origin[2])));
  // at 0 degrees
  //real p2_x = detector_x;
  //real p2_y = 0.0;
  // path length from source to detector is independent of rotation
  recon_type d_conv = recon_type(3.0 * detector_x);
  const real source_x = -2.0 * detector_x;

  // Todo - do we need this?
  recon_1d vox_z(nz + 1);
  for (int i = 0; i <= nz; i++)
    vox_z[i] = vox_origin[2] + real(i) * vox_size[2];
  real_1d yvals(ny + 1);
  for (int j = 0; j <= ny; j++)
    yvals[j] = vox_origin[1] + real(j) * vox_size[1];

  std::vector<int> mapping(nv_pixels);
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
#ifdef TEST2D
    real ymax;
    real ytop;
    if (cphi > 0.0) {
      if (sphi > 0.0) {
	ymax = y01;
	if (y11 < y00) {
	  ytop = y00;
	} else {
	  ytop = y11;
	}
      } else {
	ymax = y11;
	if (y01 < y10) {
	  ytop = y10;
	} else {
	  ytop = y01;
	}
      }
    } else {
      if (sphi > 0.0) {
	ymax = y00;
	if (y01 < y10) {
	  ytop = y10;
	} else {
	  ytop = y01;
	}
      } else {
	ymax = y10;
	if (y11 < y00) {
	  ytop = y00;
	} else {
	  ytop = y11;
	}
      }
    }
    real yx = ymax - ytop;
    if (std::abs(y_offset[a] - yx) > epsilon)
      std::cerr << "Offset " << a << ' ' << yx << ' ' << y_offset[a] << '\n';
#endif // TEST2D
    i_offset[a] = 1.0 / y_offset[a];
    if (std::abs(cphi) > std::abs(sphi))
      length[a] = vox_size[0] / std::abs(cphi);
    else
      length[a] = vox_size[1] / std::abs(sphi);
  }
  const real ihp_step = 1.0 / (h_pixels[1] - h_pixels[0]);
  const real h_pix0 = h_pixels[0] / (h_pixels[1] - h_pixels[0]);

#pragma omp parallel for shared(h_pixels, v_pixels, pixels, voxels, angles, d_conv, vox_size, vox_origin, yvals, mapping) firstprivate(n_angles, nh_pixels, nv_pixels, nx, ny, nz, detector_x, source_x, y_offset, i_offset, length, h_pix0, ihp_step, c_angle, s_angle, map_type) schedule(dynamic)
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      const real x_0 = vox_origin[0] + real(i) * vox_size[0];
      const real x_n = vox_origin[0] + real(i + 1) * vox_size[0];
      bproject_ah(source_x, detector_x, pixels, voxels,
		  x_0, yvals[j], x_n, yvals[j + 1], vox_origin[2],
		  vox_size[0], vox_size[1], vox_size[2], nx, ny, nz, i, j,
		  n_angles, nh_pixels, nv_pixels, h_pixels, v_pixels,
		  c_angle, s_angle, y_offset, i_offset, length,
		  h_pix0, ihp_step, d_conv, mapping, map_type);
    }
  }      
}
