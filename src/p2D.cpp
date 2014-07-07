
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
				    const std::vector<int> &i,
				    const std::vector<int> &j, const int n,
				    const int a, const int h,
				    const int nv, const int nz,
				    const std::vector<int> &mapping,
				    const int map_type, recon_2d &zpix)
{
  switch (map_type) {
  case 1:
    for (int m = 0; m < n; m++) {
      for (int v = 0; v < nv; v++)
	pixels[a][h][v] += voxels[i[m]][j[m]][v] * l_xy[m];
    }
#ifdef TEST3D
    for (int m = 0; m < n; m++) {
      for (int v = 0; v < nv; v++)
	zpix[m][v] += l_xy[m];
    }
#endif // TEST3D
    break;
  case 2:
    {
      int v2 = nv / 2;
      for (int m = 0; m < n; m++) {
	int v = 0;
	for (int l = 0; l < v2; l++) {
	  pixels[a][h][v] += voxels[i[m]][j[m]][l] * l_xy[m];
	  pixels[a][h][v + 1] += voxels[i[m]][j[m]][l] * l_xy[m];
	  v += 2;
	}
      }
    }
#ifdef TEST3D
    {
      int v2 = nv / 2;
      for (int m = 0; m < n; m++) {
	int v = 0;
	for (int l = 0; l < v2; l++) {
	  zpix[m][v] += l_xy[m];
	  zpix[m][v + 1] += l_xy[m];
	  v += 2;
	}
      }
    }
#endif // TEST3D
    break;
  case 4:
    {
      int v4 = nv / 4;
      for (int m = 0; m < n; m++) {
	int v = 0;
	for (int l = 0; l < v4; l++) {
	  pixels[a][h][v] += voxels[i[m]][j[m]][l] * l_xy[m];
	  pixels[a][h][v + 1] += voxels[i[m]][j[m]][l] * l_xy[m];
	  pixels[a][h][v + 2] += voxels[i[m]][j[m]][l] * l_xy[m];
	  pixels[a][h][v + 3] += voxels[i[m]][j[m]][l] * l_xy[m];
	  v += 4;
	}
      }
    }
#ifdef TEST3D
    {
      int v4 = nv / 4;
      for (int m = 0; m < n; m++) {
	int v = 0;
	for (int l = 0; l < v4; l++) {
	  zpix[m][v] += l_xy[m];
	  zpix[m][v + 1] += l_xy[m];
	  zpix[m][v + 2] += l_xy[m];
	  zpix[m][v + 3] += l_xy[m];
	  v += 4;
	}
      }
    }
#endif // TEST3D
    break;
  default:
    for (int m = 0; m < n; m++) {
      for (int v = 0; v < nv; v++)
	pixels[a][h][v] += voxels[i[m]][j[m]][mapping[v]] * l_xy[m];
    }
#ifdef TEST3D
    for (int m = 0; m < n; m++) {
      for (int v = 0; v < nv; v++)
	zpix[m][v] += l_xy[m];
    }
#endif // TEST3D
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
  std::vector<int> i_arr(2 * max_n + 1);
  std::vector<int> j_arr(2 * max_n + 1);
  int count = 0;
  if (std::abs(cphi) < epsilon) {
    if (std::abs(sphi) < epsilon) {
      // Its not a line - shouldn't happen
      //return;
    } else {
      // line parallel to y
      int i = int(std::floor((p2_x - b_x) / d_x)); // == q2_x
      if (i >= 0 and i < nx) {
	if (sphi < 0.0) {
	  for (int j = ny - 1; j >= 0; j--) {
	    l_xy[count] = d_y;
	    i_arr[count] = i;
	    j_arr[count] = j;
	    count++;
	  }
	} else {
	  for (int j = 0; j < ny; j++) {
	    l_xy[count] = d_y;
	    i_arr[count] = i;
	    j_arr[count] = j;
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
	for (int i = nx - 1; i >= 0; i--) {
	  l_xy[count] = d_x;
	  i_arr[count] = i;
	  j_arr[count] = j;
	  count++;
	}
      } else {
	for (int i = 0; i < nx; i++) {
	  l_xy[count] = d_x;
	  i_arr[count] = i;
	  j_arr[count] = j;
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
	  while (x < nx and y < ny) {
	    i_arr[count] = x;
	    j_arr[count] = y;
	    if (alpha_x[x + 1] < alpha_y[y + 1] - epsilon) {
	      l_xy[count] = (alpha_x[x + 1] - alpha_p);
	      alpha_p = alpha_x[x + 1];
	      x++;
	    } else if (alpha_y[y + 1] < alpha_x[x + 1] - epsilon) {
	      l_xy[count] = (alpha_y[y + 1] - alpha_p);
	      alpha_p = alpha_y[y + 1];
	      y++;
	    } else {
	      real mx = std::max(alpha_x[x + 1], alpha_y[y + 1]);
	      l_xy[count] = (mx - alpha_p);
	      x++;
	      y++;
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
	  while (x < nx and y >= 0) {
	    i_arr[count] = x;
	    j_arr[count] = y;
	    if (alpha_x[x + 1] < alpha_y[y] - epsilon) {
	      l_xy[count] = (alpha_x[x + 1] - alpha_p);
	      alpha_p = alpha_x[x + 1];
	      x++;
	    } else if (alpha_y[y] < alpha_x[x + 1] - epsilon) {
	      l_xy[count] = (alpha_y[y] - alpha_p);
	      alpha_p = alpha_y[y];
	      y--;
	    } else {
	      real mx = std::max(alpha_x[x + 1], alpha_y[y]);
	      l_xy[count] = (mx - alpha_p);
	      x++;
	      y--;
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
	  while (x >= 0 and y < ny) {
	    i_arr[count] = x;
	    j_arr[count] = y;
	    if (alpha_x[x] < alpha_y[y + 1] - epsilon) {
	      l_xy[count] = (alpha_x[x] - alpha_p);
	      alpha_p = alpha_x[x];
	      x--;
	    } else if (alpha_y[y + 1] < alpha_x[x] - epsilon) {
	      l_xy[count] = (alpha_y[y + 1] - alpha_p);
	      alpha_p = alpha_y[y + 1];
	      y++;
	    } else {
	      real mx = std::max(alpha_y[y + 1], alpha_x[x]);
	      l_xy[count] = (mx - alpha_p);
	      x--;
	      y++;
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
	  while (x >= 0 and y >= 0) {
	    i_arr[count] = x;
	    j_arr[count] = y;
	    if (alpha_x[x] < alpha_y[y] - epsilon) {
	      l_xy[count] = (alpha_x[x] - alpha_p);
	      alpha_p = alpha_x[x];
	      x--;
	    } else if (alpha_y[y] < alpha_x[x] - epsilon) {
	      l_xy[count] = (alpha_y[y] - alpha_p);
	      alpha_p = alpha_y[y];
	      y--;
	    } else {
	      real mx = std::max(alpha_x[x], alpha_y[y]);
	      l_xy[count] = (mx - alpha_p);
	      x--;
	      y--;
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
	    if (i_arr[k] == i and j_arr[k] == j) {
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
	std::cerr << "Not found " << a << ' ' << h << ' ' << i_arr[k]
		  << ' ' << j_arr[k] << ' ' << k << ' ' << count
		  << ' ' << l_xy[k] / d_conv << '\n';
    }
  }
#endif // TEST2D
  if (count > 0) {
    recon_2d zvec(boost::extents[count][nv]);
    for (int m = 0; m < count; m++)
      for (int v = 0; v < nv; v++)
	zvec[m][v] = 0.0;
    calc_xy_z(pixels, voxels, l_xy, i_arr, j_arr, count, a, h,
	      nv, nz, mapping, map_type, zvec);
#ifdef TEST3D
    real start[3];
    real end[3];
    start[0] = p1_x;
    start[1] = p1_y;
    end[0] = p2_x;
    end[1] = p2_y;
    for (int v = 0; v < nv; v++) {
      start[2] = v_pixels[v];
      end[2] = v_pixels[v];
      recon_type ln = 0.0;
      for (int i = 0; i < nx; i++) {
	real x_0 = b_x + real(i) * d_x;
	for (int j = 0; j < ny; j++) {
	  real y_0 = b_y + real(j) * d_y;
	  recon_type lnk = 0.0;
	  for (int k = 0; k < nz; k++) {
	    real z_0 = b_z + real(k) * d_z;
	    recon_type val = 0.0;
	    if (test_3D(start, end, x_0, y_0, z_0, d_x, d_y, d_z,
			1, 1, 1, val)) {
	      ln += val;
	      lnk += val;
	    }
	  }
	  int m;
	  for (m = 0; m < count; m++) {
	    if (i_arr[m] == i and j_arr[m] == j) {
	      if (std::abs(zvec[m][v] - lnk) > epsilon)
		std::cerr << "pix wrong " << a << ' ' << h << ' ' << v
			  << ' ' << lnk << ' ' << zvec[m][v] << ' '
			  << lnk-zvec[m][v] << ' ' << i << ' ' << j << '\n';
	      zvec[m][v] = -1.0;
	      break;
	    }
	  }
	  if (m == count and lnk > epsilon)
	    std::cerr << "pix missed " << a << ' ' << h << ' ' << v
		      << ' ' << lnk << ' ' << i << ' ' << j << '\n';
	}
      }
      for (int m = 0; m < count; m++) {
	if (zvec[m][v] > epsilon)
	  std::cerr << "pix extra " << i_arr[m] << ' ' << j_arr[m] 
		    << ' ' << zvec[m][v] << ' ' << a << ' ' << h
		    << ' ' << v << '\n';
      }
    }
#endif // TEST3D
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
  
  //#pragma omp parallel for shared(h_pixels, v_pixels, pixels, voxels, angles, d_conv, delta_z, inv_delz, vox_z, voxel_size, grid_offset) firstprivate(n_angles, n_h, n_v, nx_voxels, ny_voxels, nz_voxels) schedule(dynamic)
  for (int a = 0; a < n_angles; a++) {
    // rotate source and detector positions by current angle
    real cosa = std::cos(angles[a]);
    real sina = std::sin(angles[a]);

    for (int h = 0; h < nh_pixels; h++) {
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
				    const std::vector<int> &a,
				    const std::vector<int> &h,
				    const int n, const int i, const int j,
				    const int nv, const int nz,
				    const std::vector<int> &mapping,
				    const int map_type, recon_2d &zpix)
{
  switch (map_type) {
  case 1:
    for (int m = 0; m < n; m++) {
      for (int v = 0; v < nv; v++)
	voxels[i][j][v] += pixels[a[m]][h[m]][v] * l_xy[m];
    }
#ifdef TEST3D
    for (int m = 0; m < n; m++) {
      for (int v = 0; v < nv; v++)
	zpix[m][v] += l_xy[m];
    }
#endif // TEST3D
    break;
  case 2:
    {
      int v2 = nv / 2;
      for (int m = 0; m < n; m++) {
	int v = 0;
	for (int l = 0; l < v2; l++) {
	  voxels[i][j][l] += pixels[a[m]][h[m]][v] * l_xy[m];
	  voxels[i][j][l] += pixels[a[m]][h[m]][v + 1] * l_xy[m];
	  v += 2;
	}
      }
    }
#ifdef TEST3D
    {
      int v2 = nv / 2;
      for (int m = 0; m < n; m++) {
	int v = 0;
	for (int l = 0; l < v2; l++) {
	  zpix[m][l] += l_xy[m];
	  zpix[m][l] += l_xy[m];
	  v += 2;
	}
      }
    }
#endif // TEST3D
    break;
  case 4:
    {
      int v4 = nv / 4;
      for (int m = 0; m < n; m++) {
	int v = 0;
	for (int l = 0; l < v4; l++) {
	  voxels[i][j][l] += pixels[a[m]][h[m]][v] * l_xy[m];
	  voxels[i][j][l] += pixels[a[m]][h[m]][v + 1] * l_xy[m];
	  voxels[i][j][l] += pixels[a[m]][h[m]][v + 2] * l_xy[m];
	  voxels[i][j][l] += pixels[a[m]][h[m]][v + 3] * l_xy[m];
	  v += 4;
	}
      }
    }
#ifdef TEST3D
    {
      int v4 = nv / 4;
      for (int m = 0; m < n; m++) {
	int v = 0;
	for (int l = 0; l < v4; l++) {
	  zpix[m][l] += l_xy[m];
	  zpix[m][l] += l_xy[m];
	  zpix[m][l] += l_xy[m];
	  zpix[m][l] += l_xy[m];
	  v += 4;
	}
      }
    }
#endif // TEST3D
    break;
  default:
    for (int m = 0; m < n; m++) {
      for (int v = 0; v < nv; v++)
	voxels[i][j][mapping[v]] += pixels[a[m]][h[m]][v] * l_xy[m];
    }
#ifdef TEST3D
    for (int m = 0; m < n; m++) {
      for (int v = 0; v < nv; v++)
	zpix[m][mapping[v]] += l_xy[m];
    }
#endif // TEST3D
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
				      const recon_type d_conv,
				      const std::vector<int> &mapping,
				      const int map_type)
{
  // Rather than using the centre just calculate for all 4 corners,
  // generate h values and loop from smallest to largest.
  // Todo - check how to calculate this
  const int pix_per_vox = n_v / (nx - 1);
  // How big should the array be - Todo - use mapping for pix_per_vox?
  int count = 0;
  std::vector<int> a_arr(2 * pix_per_vox * n_angles);
  std::vector<int> h_arr(2 * pix_per_vox * n_angles);
  recon_1d l_xy(2 * pix_per_vox * n_angles);
  // corners (x0,y0), (x0,yn), (xn,y0), (xn,yn)
  // Todo - in parallel we can probably make a better guess at which 2 corners
  // we need for the upper and lower limits.
  real pixel_step = h_pixels[1] - h_pixels[0];
  for (int a = 0; a < n_angles; a++) {
    real cphi = cangle[a];
    real sphi = sangle[a];
    // Todo - is there a direction issue here that should swap 0/n?
    if (std::abs(cphi) < epsilon) {
      // line is parallel to y, so just find h that is between x_0/x_n
      if (sphi < 0.0) {
	int hmin = int(std::floor((x_0 - h_pixels[0]) / pixel_step));
	if (hmin < 0)
	  hmin = 0;
	// decrease x_n by tol so inside upper boundary, not on it?
	int hmax = int(std::floor((x_n - epsilon - h_pixels[0]) / pixel_step));
	if (hmax >= n_h)
	  hmax = n_h - 1;
	for (int h = hmin; h <= hmax; h++) {
	  if (h_pixels[h] >= x_0 and h_pixels[h] < x_n) {
	    l_xy[count] = d_y;
	    a_arr[count] = a;
	    h_arr[count] = h;
	    count++;
	  }
	}
      } else {
	int hmin = int(std::ceil((- x_n + epsilon - h_pixels[0])/pixel_step));
	if (hmin < 0)
	  hmin = 0;
	// decrease x_n by tol so inside upper boundary, not on it?
	int hmax = int(std::ceil((- x_0 - h_pixels[0]) / pixel_step));
	if (hmax >= n_h)
	  hmax = n_h - 1;
	for (int h = hmin; h <= hmax; h++) {
	  if (-h_pixels[h] >= x_0 and -h_pixels[h] < x_n) {
	    l_xy[count] = d_y;
	    a_arr[count] = a;
	    h_arr[count] = h;
	    count++;
	  }
	}
      }
    } else if (std::abs(sphi) < epsilon) {
      if (cphi < 0.0) {
	int hmin = int(std::ceil((- y_n + epsilon - h_pixels[0])/pixel_step));
	if (hmin < 0)
	  hmin = 0;
	// decrease y_n by tol so inside upper boundary, not on it?
	int hmax = int(std::ceil((- y_0 - h_pixels[0]) / pixel_step));
	if (hmax >= n_h)
	  hmax = n_h - 1;
	for (int h = hmin; h <= hmax; h++) {
	  if (- h_pixels[h] >= y_0 and - h_pixels[h] < y_n) {
	    l_xy[count] = d_x;
	    a_arr[count] = a;
	    h_arr[count] = h;
	    count++;
	  }
	}
      } else {
	int hmin = int(std::floor((y_0 - h_pixels[0]) / pixel_step));
	if (hmin < 0)
	  hmin = 0;
	// decrease y_n by tol so inside upper boundary, not on it?
	int hmax = int(std::floor((y_n - epsilon - h_pixels[0]) / pixel_step));
	if (hmax >= n_h)
	  hmax = n_h - 1;
	for (int h = hmin; h <= hmax; h++) {
	  if (h_pixels[h] >= y_0 and h_pixels[h] < y_n) {
	    l_xy[count] = d_x;
	    a_arr[count] = a;
	    h_arr[count] = h;
	    count++;
	  }
	}
      }
    } else {
      // In derivation let initial point for detector line be (dx, 0)
      real qx = detector_x * cphi; // - source_y * sphi;
      real qy = detector_x * sphi; // + source_y * cphi;
      // Then let the initial position of the be (vx, vy) not (px, py)
      // since we know the angle of all the lines which is (cphi, sphi)
      // the angle is delta_x/y
      real qpx_0 = qx - x_0;
      real qpy_0 = qy - y_0;
      real qpx_n = qx - x_n;
      real qpy_n = qy - y_n;
      /*
      real y00;
      real y01;
      real y10;
      real y11;
      if (std::abs(sphi) > std::abs(cphi)) {
	real tphi = cphi / sphi;
	real u = (qpy_0 + qpx_0 * tphi) / (sphi + cphi * tphi);
	real t = (-qpx_0 + u * cphi) / sphi;
	y00 = - t;
	u = (qpy_n + qpx_0 * tphi) / (sphi + cphi * tphi);
	t = (-qpx_0 + u * cphi) / sphi;
	y01 = - t;
	u = (qpy_0 + qpx_n * tphi) / (sphi + cphi * tphi);
	t = (-qpx_n + u * cphi) / sphi;
	y10 = - t;
	u = (qpy_n + qpx_n * tphi) / (sphi + cphi * tphi);
	t = (-qpx_n + u * cphi) / sphi;
	y11 = - t;

	real u = qpy_0 * sphi + qpx_0 * cphi;
	y00 = (qpx_0 - u * cphi) / sphi;
	u = qpy_n * sphi + qpx_0 * cphi;
	y01 = (qpx_0 - u * cphi) / sphi;
	u = qpy_0 * sphi + qpx_n * cphi;
	y10 = (qpx_n - u * cphi) / sphi;
	u = qpy_n * sphi + qpx_n * cphi;
	y11 = (qpx_n - u * cphi) / sphi;

	y00 = -qpy_0 * cphi + qpx_0 * sphi;
	y01 = -qpy_n * cphi + qpx_0 * sphi;
	y10 = -qpy_0 * cphi + qpx_n * sphi;
	y11 = -qpy_n * cphi + qpx_n * sphi;
      } else {
	real tphi = sphi / cphi;
	real u = (qpx_0 + qpy_0 * tphi) / (cphi + sphi * tphi);
	real t = (qpy_0 - u * sphi) / cphi;
	y00 = - t;
	u = (qpx_0 + qpy_n * tphi) / (cphi + sphi * tphi);
	t = (qpy_n - u * sphi) / cphi;
	y01 = - t;
	u = (qpx_n + qpy_0 * tphi) / (cphi + sphi * tphi);
	t = (qpy_0 - u * sphi) / cphi;
	y10 = - t;
	u = (qpx_n + qpy_n * tphi) / (cphi + sphi * tphi);
	t = (qpy_n - u * sphi) / cphi;
	y11 = - t;

	real u = qpx_0 * cphi + qpy_0 * sphi;
	y00 = (-qpy_0 + u * sphi) / cphi;
	u = qpx_0 * cphi + qpy_n * sphi;
	y01 = (-qpy_n + u * sphi) / cphi;
	u = qpx_n * cphi + qpy_0 * sphi;
	y10 = (- qpy_0 + u * sphi) / cphi;
	u = qpx_n * cphi + qpy_n * sphi;
	y11 = (-qpy_n + u * sphi) / cphi;

	y00 = qpx_0 * sphi - qpy_0 * cphi;
	y01 = qpx_0 * sphi - qpy_n * cphi;
	y10 = qpx_n * sphi - qpy_0 * cphi;
	y11 = qpx_n * sphi - qpy_n * cphi;
      }
      */
      real y00 = qpx_0 * sphi - qpy_0 * cphi;
      real y01 = qpx_0 * sphi - qpy_n * cphi;
      real y10 = qpx_n * sphi - qpy_0 * cphi;
      real y11 = qpx_n * sphi - qpy_n * cphi;
      int h00 = int(std::floor((y00 - h_pixels[0]) / pixel_step));
      int h01 = int(std::floor((y01 - h_pixels[0]) / pixel_step));
      int hmin = std::min(h00, h01);
      int hmax = std::max(h00, h01);
      int h10 = int(std::floor((y10 - h_pixels[0]) / pixel_step));
      hmin = std::min(hmin, h10);
      hmax = std::max(hmax, h10);
      int h11 = int(std::floor((y11 - h_pixels[0]) / pixel_step));
      hmin = std::max(std::min(hmin, h11), 0);
      hmax = std::min(std::max(hmax, h11), n_h - 1);
      // If it intercepts voxel then calc length and store
      /**/
      // Use an absolute paramteric form from 0 to L, not 0 to 1
      // I think the min -L, max 0 limits are correct
      const real delta_x = cphi;
      const real delta_y = sphi;
      const real inv_dx = 1.0 / delta_x;
      const real inv_dy = 1.0 / delta_y;
      for (int h = hmin; h <= hmax; h++) {
	const real p2_x = cphi * detector_x - sphi * h_pixels[h];
	const real p2_y = sphi * detector_x + cphi * h_pixels[h];
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
	  l_xy[count] = (alpha_max - alpha_min);
	  a_arr[count] = a;
	  h_arr[count] = h;
	  count++;
	}
      }
      /**/
      /*
      for (int h = hmin; h <= hmax; h++) {
	const real p2_x = cphi * detector_x - sphi * h_pixels[h];
	const real p2_y = sphi * detector_x + cphi * h_pixels[h];
	const real p1_x = p2_x - real(3.0) * cphi * detector_x;
	const real p1_y = p2_y - real(3.0) * sphi * detector_x;
	const real delta_x = p2_x - p1_x; // cphi - Todo
	const real delta_y = p2_y - p1_y; // sphi - Todo
	const real inv_dx = 1.0 / delta_x;
	const real inv_dy = 1.0 / delta_y;
	const real alpha_x_0 = (x_0 - p1_x) * inv_dx;
	const real alpha_y_0 = (y_0 - p1_y) * inv_dy;
	const real alpha_x_n = (x_n - p1_x) * inv_dx;
	const real alpha_y_n = (y_n - p1_y) * inv_dy;
	const real alpha_x_min = std::min(alpha_x_0, alpha_x_n);
	const real alpha_x_max = std::max(alpha_x_0, alpha_x_n);
	const real alpha_y_min = std::min(alpha_y_0, alpha_y_n);
	const real alpha_y_max = std::max(alpha_y_0, alpha_y_n);
	const real alpha_min = std::max(std::max(alpha_x_min, alpha_y_min),
					real(0.0));
	const real alpha_max = std::min(std::min(alpha_x_max, alpha_y_max),
					real(1.0));
	if (alpha_min < alpha_max - epsilon) {
	  l_xy[count] = (alpha_max - alpha_min) * d_conv;
	  a_arr[count] = a;
	  h_arr[count] = h;
	  count++;
	}
      }
      */
    }
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
	    if (a_arr[k] == a and h_arr[k] == h) {
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
	std::cerr << "Not found " << i << ' ' << j << ' ' << a_arr[k]
		  << ' ' << h_arr[k] << ' ' << k << ' ' << count
		  << ' ' << l_xy[k] / d_conv << '\n';
    }
  }
#endif // TEST2D
  if (count > 0) {
    recon_2d zvec(boost::extents[count][nz]);
    for (int m = 0; m < count; m++)
      for (int k = 0; k < nz; k++)
	zvec[m][k] = 0.0;
    calc_ah_z(pixels, voxels, l_xy, a_arr, h_arr, count,
	      i, j, n_v, nz, mapping, map_type, zvec);
#ifdef TEST3D
    real start[3];
    real end[3];
    for (int k = 0; k < nz; k++) {
      real z_0 = b_z + real(k) * d_z;
      recon_type ln = 0.0;
      for (int a = 0; a < n_angles; a++) {
	for (int h = 0; h < n_h; h++) {
	  real cos_curr_angle = cangle[a];
	  real sin_curr_angle = sangle[a];
	  end[0] = cos_curr_angle * detector_x - sin_curr_angle * h_pixels[h];
	  end[1] = sin_curr_angle * detector_x + cos_curr_angle * h_pixels[h];
	  start[0] = end[0] - real(3.0) * cos_curr_angle * detector_x;
	  start[1] = end[1] - real(3.0) * sin_curr_angle * detector_x;
	  recon_type lnk = 0.0;
	  for (int v = 0; v < n_v; v++) {
	    start[2] = v_pixels[v];
	    end[2] = v_pixels[v];
	    recon_type val = 0.0;
	    if (test_3D(start, end, x_0, y_0, z_0, d_x, d_y, d_z,
			1, 1, 1, val)) {
	      ln += val;
	      lnk += val;
	    }
	  }
	  int m;
	  for (m = 0; m < count; m++) {
	    if (a_arr[m] == a and h_arr[m] == h) {
	      if (std::abs(zvec[m][k] - lnk) > epsilon)
		std::cerr << "vox wrong " << a << ' ' << h
			  << ' ' << lnk << ' ' << zvec[m][k] << ' '
			  << lnk-zvec[m][k] << ' ' << i << ' ' << j
			  << ' ' << k << '\n';
	      zvec[m][k] = -1.0;
	      break;
	    }
	  }
	  if (m == count and lnk > epsilon)
	    std::cerr << "vox missed " << a << ' ' << h
		      << ' ' << lnk << ' ' << ' ' << i << ' ' << j
		      << ' ' << k << '\n';
	}
      }
      for (int m = 0; m < count; m++) {
	if (zvec[m][k] > epsilon)
	  std::cerr << "vox extra " << a_arr[m] << ' ' << h_arr[m] 
		    << ' ' << zvec[m][k] << ' ' << i << ' ' << j
		    << ' ' << k << '\n';
      }
    }
#endif // TEST3D
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

  // #pragma
  for (int i = 0; i < nx; i++) {
    const real x_0 = vox_origin[0] + real(i) * vox_size[0];
    const real x_n = vox_origin[0] + real(i + 1) * vox_size[0];
    for (int j = 0; j < ny; j++) {
      bproject_ah(source_x, detector_x, pixels, voxels,
		  x_0, yvals[j], x_n, yvals[j + 1], vox_origin[2],
		  vox_size[0], vox_size[1], vox_size[2], nx, ny, nz, i, j,
		  n_angles, nh_pixels, nv_pixels, h_pixels, v_pixels,
		  c_angle, s_angle, d_conv, mapping, map_type);
    }
  }      
}
