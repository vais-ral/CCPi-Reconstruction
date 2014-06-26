
#include <float.h>
#ifdef MATLAB_MEX_FILE
#  include "mex_types.hpp"
#else
#  include "base_types.hpp"
#endif // mex
#include "instruments.hpp"
#include "timer.hpp"
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

void CCPi::cone_beam::calc_xy_z(pixel_data &pixels, voxel_data &voxels,
				const recon_1d &alpha_xy,
				const std::vector<int> &i,
				const std::vector<int> &j, const int n,
				const int a, const int h, const recon_type p1_z,
				const recon_type b_z, const recon_type d_z,
				const int nv, const int nz, const int midp,
				const recon_2d &d_conv, const recon_1d &delta_z,
				const recon_1d &inv_delz, const recon_1d &vox_z,
				recon_2d &zpix)
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
  const int nzm1 = nz - 1;
  const recon_type inv_dz = 1.0 / d_z;
  const recon_type pzbz = (p1_z - b_z) * inv_dz;
  for (int m = 1; m < n; m++) {
    const recon_type alpha_inv = alpha_xy[m - 1] * inv_dz;
    for (int v = midp - 1; v >= 0; v--) {
      int k = int(std::floor(pzbz + alpha_inv * delta_z[v]));
      //k = k - 1;
      if (k > 0) {
	recon_type alpha_z = vox_z[k] * inv_delz[v];
	recon_type min_z = std::min(alpha_z, alpha_xy[m]);
	pixels[a][h][v] += voxels[i[m]][j[m]][k] * (min_z - alpha_xy[m - 1]);
	pixels[a][h][v] += voxels[i[m]][j[m]][k - 1] * (alpha_xy[m] - min_z);
#ifdef TEST3D
	zpix[m][v] += (min_z - alpha_xy[m - 1]) + (alpha_xy[m] - min_z);
#endif // TEST3D
      } else if (k == 0) {
	recon_type alpha_z = vox_z[k] * inv_delz[v];
	recon_type min_z = std::min(alpha_z, alpha_xy[m]);
	pixels[a][h][v] += voxels[i[m]][j[m]][k] * (min_z - alpha_xy[m - 1]);
#ifdef TEST3D
	zpix[m][v] += (min_z - alpha_xy[m - 1]);
#endif // TEST3D
      } else
	break;
    }
    for (int v = midp; v < nv; v++) {
      int k = int(std::floor(pzbz + alpha_inv * delta_z[v]));
      if (k < nzm1) {
	recon_type alpha_z = vox_z[k + 1] * inv_delz[v];
	recon_type min_z = std::min(alpha_z, alpha_xy[m]);
	pixels[a][h][v] += voxels[i[m]][j[m]][k] * (min_z - alpha_xy[m - 1]);
	pixels[a][h][v] += voxels[i[m]][j[m]][k + 1] * (alpha_xy[m] - min_z);
#ifdef TEST3D
	zpix[m][v] += (min_z - alpha_xy[m - 1]) + (alpha_xy[m] - min_z);
#endif // TEST3D
      } else if (k == nzm1) {
	recon_type alpha_z = vox_z[k + 1] * inv_delz[v];
	recon_type min_z = std::min(alpha_z, alpha_xy[m]);
	pixels[a][h][v] += voxels[i[m]][j[m]][k] * (min_z - alpha_xy[m - 1]);
#ifdef TEST3D
	zpix[m][v] += (min_z - alpha_xy[m - 1]);
#endif // TEST3D
      } else
	break;
    }
  }
  for (int v = 0; v < nv; v++)
    pixels[a][h][v] *= d_conv[h][v];
#ifdef TEST3D
  for (int m = 1; m < n; m++)
    for (int v = 0; v < nv; v++)
      zpix[m][v] *= d_conv[h][v];
#endif // TEST3D
}

void CCPi::cone_beam::fproject_xy(const real p1_x, const real p1_y,
				  const real p2_x, const real p2_y,
				  pixel_data &pixels, voxel_data &voxels,
				  const real b_x, const real b_y,
				  const real b_z, const real d_x,
				  const real d_y, const real d_z,
				  const int nx, const int ny, const int nz,
				  const int a, const int h,
				  const real source_z, const int nv,
				  const int midp, const recon_2d &d_conv,
				  const recon_1d &delta_z,
				  const recon_1d &inv_delz,
				  const recon_1d &vox_z,
				  const real_1d &v_pixels)
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

  int max_n = std::max(nx, ny);
  recon_1d alpha_xy(2 * max_n);
  std::vector<int> i_arr(2 * max_n + 1);
  std::vector<int> j_arr(2 * max_n + 1);
  int count = 0;
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
	  i_arr[count] = i;
	  j_arr[count] = ny - 1;
	  count++;
	  for (int j = ny - 1; j >= 0; j--) {
	    alpha_xy[count] = (real(j) - q1_y) * inv_dy;
	    i_arr[count] = i;
	    j_arr[count] = j;
	    count++;
	  }
	} else {
	  alpha_xy[count] = (real(0) - q1_y) * inv_dy;
	  i_arr[count] = i;
	  j_arr[count] = 0;
	  count++;
	  for (int j = 0; j < ny; j++) {
	    alpha_xy[count] = (real(j + 1) - q1_y) * inv_dy;
	    i_arr[count] = i;
	    j_arr[count] = j;
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
	i_arr[count] = nx - 1;
	j_arr[count] = j;
	count++;
	for (int i = nx - 1; i >= 0; i--) {
	  alpha_xy[count] = (real(i) - q1_x) * inv_dx;
	  i_arr[count] = i;
	  j_arr[count] = j;
	  count++;
	}
      } else {
	alpha_xy[count] = (real(0) - q1_x) * inv_dx;
	i_arr[count] = 0;
	j_arr[count] = j;
	count++;
	for (int i = 0; i < nx; i++) {
	  alpha_xy[count] = (real(i + 1) - q1_x) * inv_dx;
	  i_arr[count] = i;
	  j_arr[count] = j;
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
	  i_arr[count] = x;
	  j_arr[count] = y;
	  count++;
	  // could do x_next/y_next here and only calc the one that changes
	  // inside the if statements below, which would reduce the flops
	  while (x < nx and y < ny) {
	    // could calc a general alpha_step and maybe have rounding issues
	    i_arr[count] = x;
	    j_arr[count] = y;
	    if (alpha_x[x + 1] < alpha_y[y + 1] - epsilon) {
	      alpha_xy[count] = alpha_x[x + 1];
	      x++;
	    } else if (alpha_y[y + 1] < alpha_x[x + 1] - epsilon) {
	      alpha_xy[count] = alpha_y[y + 1];
	      y++;
	    } else {
	      alpha_xy[count] = std::max(alpha_x[x + 1], alpha_y[y + 1]);
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
	  i_arr[count] = x;
	  j_arr[count] = y;
	  count++;
	  while (x < nx and y >= 0) {
	    i_arr[count] = x;
	    j_arr[count] = y;
	    if (alpha_x[x + 1] < alpha_y[y] - epsilon) {
	      alpha_xy[count] = alpha_x[x + 1];
	      x++;
	    } else if (alpha_y[y] < alpha_x[x + 1] - epsilon) {
	      alpha_xy[count] = alpha_y[y];
	      y--;
	    } else {
	      alpha_xy[count] = std::max(alpha_x[x + 1], alpha_y[y]);
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
	  i_arr[count] = x;
	  j_arr[count] = y;
	  count++;
	  while (x >= 0 and y < ny) {
	    i_arr[count] = x;
	    j_arr[count] = y;
	    if (alpha_x[x] < alpha_y[y + 1] - epsilon) {
	      alpha_xy[count] = alpha_x[x];
	      x--;
	    } else if (alpha_y[y + 1] < alpha_x[x] - epsilon) {
	      alpha_xy[count] = alpha_y[y + 1];
	      y++;
	    } else {
	      alpha_xy[count] = std::max(alpha_y[y + 1], alpha_x[x]);
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
	  i_arr[count] = x;
	  j_arr[count] = y;
	  count++;
	  while (x >= 0 and y >= 0) {
	    i_arr[count] = x;
	    j_arr[count] = y;
	    if (alpha_x[x] < alpha_y[y] - epsilon) {
	      alpha_xy[count] = alpha_x[x];
	      x--;
	    } else if (alpha_y[y] < alpha_x[x] - epsilon) {
	      alpha_xy[count] = alpha_y[y];
	      y--;
	    } else {
	      alpha_xy[count] = std::max(alpha_x[x], alpha_y[y]);
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
	  for (k = 1; k < count; k++) {
	    if (i_arr[k] == i and j_arr[k] == j) {
	      cmp[k] = true;
	      real diff = alpha_xy[k] - alpha_xy[k-1];
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
    for (int k = 1; k < count; k++) {
      if (!cmp[k])
	std::cerr << "Not found " << a << ' ' << h << ' ' << i_arr[k]
		  << ' ' << j_arr[k] << ' ' << k << ' ' << count
		  << ' ' << alpha_xy[k] - alpha_xy[k-1] << '\n';
    }
  }
#endif // TEST2D
  if (count > 0) {
    recon_2d zvec(boost::extents[count][nv]);
    for (int m = 0; m < count; m++)
      for (int v = 0; v < nv; v++)
	zvec[m][v] = 0.0;
    calc_xy_z(pixels, voxels, alpha_xy, i_arr, j_arr, count, a, h, source_z,
	      b_z, d_z, nv, nz, midp, d_conv, delta_z, inv_delz, vox_z, zvec);
#ifdef TEST3D
    real start[3];
    real end[3];
    start[0] = p1_x;
    start[1] = p1_y;
    start[2] = source_z;
    end[0] = p2_x;
    end[1] = p2_y;
    for (int v = 0; v < nv; v++) {
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
	  for (m = 1; m < count; m++) {
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
      for (int m = 1; m < count; m++) {
	if (zvec[m][v] > epsilon)
	  std::cerr << "pix extra " << i_arr[m] << ' ' << j_arr[m] 
		    << ' ' << zvec[m][v] << ' ' << a << ' ' << h
		    << ' ' << v << '\n';
      }
      /*
      if (std::abs(ln) > epsilon) {
	if (std::abs(ln - zvec[v]) / ln > 10.0 * epsilon)
	  std::cerr << "pix wrong " << a << ' ' << h << ' ' << v
		    << ' ' << ln << ' ' << zvec[v] << ' ' << ln-zvec[v] << '\n';
      } else if (std::abs(zvec[v]) > epsilon) {
	if (std::abs(ln - zvec[v]) / zvec[v] > 10.0 * epsilon)
	  std::cerr << "pix wrong " << a << ' ' << h << ' ' << v
		    << ' ' << ln << ' ' << zvec[v] << ' ' << ln-zvec[v] << '\n';
      }
      */
    }
#endif // TEST3D
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

  recon_1d delta_z(n_v);
  for (int i = 0; i < n_v; i++)
    delta_z[i] = v_pixels[i] - source_z;
  recon_1d inv_delz(n_v);
  for (int i = 0; i < n_v; i++)
    inv_delz[i] = 1.0 / delta_z[i];
  recon_1d vox_z(nz_voxels + 1);
  for (int i = 0; i <= nz_voxels; i++)
    vox_z[i] = grid_offset[2] + real(i) * voxel_size[2] - source_z;
  // Todo - 1d arrays of x,y positions for voxels?


  //#pragma omp parallel for shared(h_pixels, v_pixels, pixels, voxels, angles, d_conv, delta_z, inv_delz, vox_z, voxel_size, grid_offset) firstprivate(n_angles, n_h, n_v, nx_voxels, ny_voxels, nz_voxels) schedule(dynamic)
  for (int a = 0; a < n_angles; a++) {
    /* rotate source and detector positions by current angle */
    real cosa = std::cos(angles[a]);
    real sina = std::sin(angles[a]);
    real p1_x = cosa * source_x - sina * source_y;
    real p1_y = sina * source_x + cosa * source_y;

    for (int h = 0; h < n_h; h++) {
      real p2_x = cosa * detector_x - sina * h_pixels[h];
      real p2_y = sina * detector_x + cosa * h_pixels[h];
      fproject_xy(p1_x, p1_y, p2_x, p2_y, pixels, voxels, grid_offset[0],
		  grid_offset[1], grid_offset[2], voxel_size[0], voxel_size[1],
		  voxel_size[2], nx_voxels, ny_voxels, nz_voxels, a, h,
		  source_z, n_v, mid, d_conv, delta_z, inv_delz, vox_z,
		  v_pixels);
    }
  }
}

// Todo - loop over k and find v? alpha_xy may be the wrong info for this.
void CCPi::cone_beam::calc_ah_z(pixel_data &pixels, voxel_data &voxels,
				const recon_1d &alpha_xy_0,
				const recon_1d &alpha_xy_1,
				const std::vector<int> &a,
				const std::vector<int> &h,
				const int n, const int i, const int j,
				const recon_type p1_z, const recon_type b_z,
				const recon_type d_z, const int nv,
				const int nz, const int midp,
				const recon_2d &d_conv, const recon_1d &delta_z,
				const recon_1d &inv_delz, const recon_1d &vox_z,
				recon_2d &zpix)
{
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
  const int nzm1 = nz - 1;
  const recon_type inv_dz = 1.0 / d_z;
  const recon_type pzbz = (p1_z - b_z) * inv_dz;
  for (int m = 0; m < n; m++) {
    const recon_type alpha_inv = alpha_xy_0[m] * inv_dz;
    for (int v = midp - 1; v >= 0; v--) {
      // truncation of -0.2 -> 0 issue, when we want to terminate as k < 0
      int k = int(std::floor(pzbz + alpha_inv * delta_z[v]));
      k = k - 1;
      if (k > 0) {
	recon_type alpha_z = vox_z[k] * inv_delz[v];
	recon_type min_z = std::min(alpha_z, alpha_xy_1[m]);
	voxels[i][j][k] += pixels[a[m]][h[m]][v] * (min_z - alpha_xy_0[m])
	  * d_conv[h[m]][v];
	voxels[i][j][k - 1] += pixels[a[m]][h[m]][v] * (alpha_xy_1[m] - min_z)
	  * d_conv[h[m]][v];
#ifdef TEST3D
	zpix[m][k] += (min_z - alpha_xy_0[m]) * d_conv[h[m]][v];
	zpix[m][k - 1] += (alpha_xy_1[m] - min_z) * d_conv[h[m]][v];
#endif // TEST3D
      } else if (k == 0) {
	recon_type alpha_z = vox_z[k] * inv_delz[v];
	recon_type min_z = std::min(alpha_z, alpha_xy_1[m]);
	voxels[i][j][k] += pixels[a[m]][h[m]][v] * (min_z - alpha_xy_0[m])
	  * d_conv[h[m]][v];
#ifdef TEST3D
	zpix[m][k] += (min_z - alpha_xy_0[m]) * d_conv[h[m]][v];
#endif // TEST3D
      } else
	break;
    }
    for (int v = midp; v < nv; v++) {
      int k = int(std::floor(pzbz + alpha_inv * delta_z[v]));
      if (k < nzm1) {
	recon_type alpha_z = vox_z[k + 1] * inv_delz[v];
	recon_type min_z = std::min(alpha_z, alpha_xy_1[m]);
	voxels[i][j][k] += pixels[a[m]][h[m]][v] * (min_z - alpha_xy_0[m])
	  * d_conv[h[m]][v];
	voxels[i][j][k + 1] += pixels[a[m]][h[m]][v] * (alpha_xy_1[m] - min_z)
	  * d_conv[h[m]][v];
#ifdef TEST3D
	zpix[m][k] += (min_z - alpha_xy_0[m]) * d_conv[h[m]][v];
	zpix[m][k + 1] += (alpha_xy_1[m] - min_z) * d_conv[h[m]][v];
#endif // TEST3D
      } else if (k == nzm1) {
	recon_type alpha_z = vox_z[k + 1] * inv_delz[v];
	recon_type min_z = std::min(alpha_z, alpha_xy_1[m]);
	voxels[i][j][k] += pixels[a[m]][h[m]][v] * (min_z - alpha_xy_0[m])
	  * d_conv[h[m]][v];
#ifdef TEST3D
	zpix[m][k] += (min_z - alpha_xy_0[m]) * d_conv[h[m]][v];
#endif // TEST3D
      } else
	break;
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
				  const real detector_x, pixel_data &pixels,
				  voxel_data &voxels, const real b_x,
				  const real b_y, const real b_z,
				  const real d_x, const real d_y,
				  const real d_z, const int nx, const int ny,
				  const int nz, const int i, const int j,
				  const real source_z, const int n_angles,
				  const int n_h, const int n_v,
				  const real_1d &h_pixels,
				  const real_1d &v_pixels, const int midp,
				  const real_1d &cangle, const real_1d &sangle,
				  const recon_2d &d_conv,
				  const recon_1d &delta_z,
				  const recon_1d &inv_delz,
				  const recon_1d &vox_z)
{
  // Rather than using the centre just calculate for all 4 corners,
  // generate h values and loop from smallest to largest.
  // Todo - check how to calculate this
  const int pix_per_vox = n_v / (nx - 1);
  // How big should the array be
  int count = 0;
  std::vector<int> a_arr(2 * pix_per_vox * n_angles);
  std::vector<int> h_arr(2 * pix_per_vox * n_angles);
  recon_1d alpha_xy_0(2 * pix_per_vox * n_angles);
  recon_1d alpha_xy_1(2 * pix_per_vox * n_angles);
  // corners (x0,y0), (x0,yn), (xn,y0), (xn,yn)
  const real x_0 = b_x + real(i) * d_x;
  const real y_0 = b_y + real(j) * d_y;
  const real x_n = b_x + real(i + 1) * d_x;
  const real y_n = b_y + real(j + 1) * d_y;
  real pixel_step = h_pixels[1] - h_pixels[0];
  for (int a = 0; a < n_angles; a++) {
    real cphi = cangle[a];
    real sphi = sangle[a];
    real px = source_x * cphi - source_y * sphi;
    real py = source_x * sphi + source_y * cphi;
    real qx = detector_x * cphi - source_y * sphi;
    real qy = detector_x * sphi + source_y * cphi;
    real delta_x = x_0 - px;
    real delta_y = y_0 - py;
    real t;
    if (std::abs(sphi) > std::abs(cphi)) {
      real u = (qy - py + (qx - px) * cphi / sphi) / 
	(delta_y + delta_x * cphi / sphi);
      t = (px - qx + u * delta_x) / sphi;
    } else {
      real u = (qx - px + (qy - py) * sphi / cphi) /
	(delta_x + delta_y * sphi / cphi);
      t = (qy - py - u * delta_y) / cphi;
    }
    real y00 = source_y - t;
    int h00 = int(std::floor((y00 - h_pixels[0]) / pixel_step));
    delta_x = x_0 - px;
    delta_y = y_n - py;
    if (std::abs(sphi) > std::abs(cphi)) {
      real u = (qy - py + (qx - px) * cphi / sphi) / 
	(delta_y + delta_x * cphi / sphi);
      t = (px - qx + u * delta_x) / sphi;
    } else {
      real u = (qx - px + (qy - py) * sphi / cphi) /
	(delta_x + delta_y * sphi / cphi);
      t = (qy - py - u * delta_y) / cphi;
    }
    real y01 = source_y - t;
    int h01 = int(std::floor((y01 - h_pixels[0]) / pixel_step));
    int hmin = std::min(h00, h01);
    int hmax = std::max(h00, h01);
    delta_x = x_n - px;
    delta_y = y_0 - py;
    if (std::abs(sphi) > std::abs(cphi)) {
      real u = (qy - py + (qx - px) * cphi / sphi) / 
	(delta_y + delta_x * cphi / sphi);
      t = (px - qx + u * delta_x) / sphi;
    } else {
      real u = (qx - px + (qy - py) * sphi / cphi) /
	(delta_x + delta_y * sphi / cphi);
      t = (qy - py - u * delta_y) / cphi;
    }
    real y10 = source_y - t;
    int h10 = int(std::floor((y10 - h_pixels[0]) / pixel_step));
    hmin = std::min(hmin, h10);
    hmax = std::max(hmax, h10);
    delta_x = x_n - px;
    delta_y = y_n - py;
    if (std::abs(sphi) > std::abs(cphi)) {
      real u = (qy - py + (qx - px) * cphi / sphi) / 
	(delta_y + delta_x * cphi / sphi);
      t = (px - qx + u * delta_x) / sphi;
    } else {
      real u = (qx - px + (qy - py) * sphi / cphi) /
	(delta_x + delta_y * sphi / cphi);
      t = (qy - py - u * delta_y) / cphi;
    }
    real y11 = source_y - t;
    int h11 = int(std::floor((y11 - h_pixels[0]) / pixel_step));
    hmin = std::max(std::min(hmin, h11), 0);
    hmax = std::min(std::max(hmax, h11), n_h - 1);
    // If it intercepts voxel then calc alpha_xy_0, alpha_xy_1 and store
    const real p1_x = px;
    const real p1_y = py;
    for (int h = hmin; h <= hmax; h++) {
      const real p2_x = cphi * detector_x - sphi * h_pixels[h];
      const real p2_y = sphi * detector_x + cphi * h_pixels[h];
      const real delta_x = p2_x - p1_x;
      const real delta_y = p2_y - p1_y;
      if (std::abs(delta_x) < epsilon) {
	report_error("Ooops - delta x");
      } else if (std::abs(delta_y) < epsilon) {
	report_error("Ooops - delta y");
      } else {
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
	  alpha_xy_0[count] = alpha_min;
	  alpha_xy_1[count] = alpha_max;
	  a_arr[count] = a;
	  h_arr[count] = h;
	  count++;
	}
      }
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
	start[0] = cos_curr_angle * source_x - sin_curr_angle * source_y;
	start[1] = sin_curr_angle * source_x + cos_curr_angle * source_y;
	end[0] = cos_curr_angle * detector_x - sin_curr_angle * h_pixels[h];
	end[1] = sin_curr_angle * detector_x + cos_curr_angle * h_pixels[h];
	recon_type ln1, ln2;
	if (test_2D(start, end, x_0, y_0, d_x, d_y, 1, 1, ln1, ln2)) {
	  recon_type ln = ln1 - ln2;
	  int k;
	  for (k = 0; k < count; k++) {
	    if (a_arr[k] == a and h_arr[k] == h) {
	      cmp[k] = true;
	      real diff = alpha_xy_1[k] - alpha_xy_0[k];
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
		  << ' ' << alpha_xy_1[k] - alpha_xy_0[k] << '\n';
    }
  }
#endif // TEST2D
  if (count > 0) {
    recon_2d zvec(boost::extents[count][nz]);
    for (int m = 0; m < count; m++)
      for (int k = 0; k < nz; k++)
	zvec[m][k] = 0.0;
    calc_ah_z(pixels, voxels, alpha_xy_0, alpha_xy_1, a_arr, h_arr, count,
	      i, j, source_z, b_z, d_z, n_v, nz, midp, d_conv, delta_z,
	      inv_delz, vox_z, zvec);
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
	  start[0] = cos_curr_angle * source_x - sin_curr_angle * source_y;
	  start[1] = sin_curr_angle * source_x + cos_curr_angle * source_y;
	  start[2] = source_z;
	  /* loop over y values on detector */
	  end[0] = cos_curr_angle * detector_x - sin_curr_angle * h_pixels[h];
	  end[1] = sin_curr_angle * detector_x + cos_curr_angle * h_pixels[h];
	  recon_type lnk = 0.0;
	  for (int v = 0; v < n_v; v++) {
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
      /*
      if (std::abs(ln) > epsilon) {
	if (std::abs(ln - zvec[k]) / ln > 10.0 * epsilon)
	  std::cerr << "vox wrong " << i << ' ' << j << ' ' << k
		    << ' ' << ln << ' ' << zvec[k] << ' ' << ln-zvec[k] << '\n';
      } else if (std::abs(zvec[k]) > epsilon) {
	if (std::abs(ln - zvec[k]) / zvec[k] > 10.0 * epsilon)
	  std::cerr << "vox wrong " << i << ' ' << j << ' ' << k
		    << ' ' << ln << ' ' << zvec[k] << ' ' << ln-zvec[k] << '\n';
      }
      */
    }
#endif // TEST3D
  }
}

void CCPi::cone_beam::b2D(const real source_x, const real source_y,
			  const real source_z, const real detector_x,
			  const real_1d &h_pixels, const real_1d &v_pixels,
			  const real_1d &angles, pixel_data &pixels,
			  voxel_data &voxels, const int n_angles,
			  const int n_h, const int n_v,
			  const real vox_origin[3], const real vox_size[3],
			  const int nx, const int ny, const int nz)
{
  // Todo - check that source_z to det_z max z angle < 45 for safe usage.
  int mid = -1;
  for (int i = 0; i < n_v; i++) {
    if (v_pixels[i] > source_z) {
      mid = i;
      break;
    }
  }

  // Todo - 1d arrays of x,y,p1x.p1y positions and 2d p2x/p2y?

  //recon_2d line_angles(boost::extents[n_angles][n_h]);
  real_1d c_angle(n_angles);
  real_1d s_angle(n_angles);
  for (int a = 0; a < n_angles; a++) {
    real cos_phi = std::cos(angles[a]);
    real sin_phi = std::sin(angles[a]);
    c_angle[a] = cos_phi;
    s_angle[a] = sin_phi;
  }
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

  recon_1d delta_z(n_v);
  for (int i = 0; i < n_v; i++)
    delta_z[i] = v_pixels[i] - source_z;
  recon_1d inv_delz(n_v);
  for (int i = 0; i < n_v; i++)
    inv_delz[i] = 1.0 / delta_z[i];
  recon_1d vox_z(nz + 1);
  for (int i = 0; i <= nz; i++)
    vox_z[i] = vox_origin[2] + real(i) * vox_size[2] - source_z;

  // #pragma
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      bproject_ah(source_x, source_y, detector_x, pixels, voxels,
		  vox_origin[0], vox_origin[1], vox_origin[2],
		  vox_size[0], vox_size[1], vox_size[2], nx, ny, nz, i, j,
		  source_z, n_angles, n_h, n_v, h_pixels, v_pixels, mid,
		  c_angle, s_angle, d_conv, delta_z, inv_delz, vox_z);
    }
  }      
}
