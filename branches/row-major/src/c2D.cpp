
#ifdef MATLAB_MEX_FILE
#  include "mex_types.hpp"
#else
#  include "base_types.hpp"
#endif // mex
#include "instruments.hpp"
#include "timer.hpp"
#include "ui_calls.hpp"

void CCPi::cone_beam::calc_xy_z(pixel_data &pixels, voxel_data &voxels,
				const recon_1d &alpha_xy,
				const std::vector<int> &i,
				const std::vector<int> &j, const int n,
				const int a, const int h, const recon_type p1_z,
				const recon_type b_z, const recon_type d_z,
				const int nv, const int nz, const int midp,
				const recon_2d &d_conv, const recon_1d &delta_z,
				const recon_1d &inv_delz, const recon_1d &vox_z)
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
      int k = int(pzbz + alpha_inv * delta_z[v]);
      if (k > 0) {
	recon_type alpha_z = vox_z[k] * inv_delz[v];
	recon_type min_z = std::min(alpha_z, alpha_xy[m]);
	pixels[a][h][v] += voxels[i[m]][j[m]][k] * (min_z - alpha_xy[m - 1]);
	pixels[a][h][v] += voxels[i[m]][j[m]][k - 1] * (alpha_xy[m] - min_z);
      } else if (k == 0) {
	recon_type alpha_z = vox_z[k] * inv_delz[v];
	recon_type min_z = std::min(alpha_z, alpha_xy[m]);
	pixels[a][h][v] += voxels[i[m]][j[m]][k] * (min_z - alpha_xy[m - 1]);
      } else
	break;
    }
    for (int v = midp; v < nv; v++) {
      int k = int(pzbz + alpha_inv * delta_z[v]);
      if (k < nzm1) {
	recon_type alpha_z = vox_z[k + 1] * inv_delz[v];
	recon_type min_z = std::min(alpha_z, alpha_xy[m]);
	pixels[a][h][v] += voxels[i[m]][j[m]][k] * (min_z - alpha_xy[m - 1]);
	pixels[a][h][v] += voxels[i[m]][j[m]][k + 1] * (alpha_xy[m] - min_z);
      } else if (k == nzm1) {
	recon_type alpha_z = vox_z[k + 1] * inv_delz[v];
	recon_type min_z = std::min(alpha_z, alpha_xy[m]);
	pixels[a][h][v] += voxels[i[m]][j[m]][k] * (min_z - alpha_xy[m - 1]);
      } else
	break;
    }
  }
  for (int v = 0; v < nv; v++)
    pixels[a][h][v] *= d_conv[h][v];
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
				  const recon_1d &vox_z)
{
  const real tol = 1e-8; // use C++ limits FLT/DBL_EPSILON float.h
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
  if (std::abs(delta_x) < tol) {
    if (std::abs(delta_y) < tol) {
      // Its not a line - shouldn't happen
      //return;
    } else {
      // line parallel to y
      int i = int(q1_x); // == q2_x
      if (i >= 0 and i < nx) {
	if (delta_y < 0.0) {
	  for (int j = ny; j >= 0; j--) {
	    alpha_xy[count] = (real(j) - q1_y) * inv_dy;
	    i_arr[j] = i;
	    j_arr[j] = j;
	    count++;
	  }
	} else {
	  for (int j = 0; j <= ny; j++) {
	    alpha_xy[count] = (real(j) - q1_y) * inv_dy;
	    i_arr[j] = i;
	    j_arr[j] = j;
	    count++;
	  }
	}
      }
    }
  } else if (std::abs(delta_y) < tol) {
    // line parallel to x
    int j = int(q1_y); // == q2_y
    if (j >= 0 and j < ny) {
      if (delta_x < 0.0) {
	for (int i = nx; i >= 0; i--) {
	  alpha_xy[count] = (real(i) - q1_x) * inv_dx;
	  i_arr[i] = i;
	  j_arr[i] = j;
	  count++;
	}
      } else {
	for (int i = 0; i <= nx; i++) {
	  alpha_xy[count] = (real(i) - q1_x) * inv_dx;
	  i_arr[i] = i;
	  j_arr[i] = j;
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
    
    if (alpha_min < alpha_max) {
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
	    y = int(q1_y + alpha_min * delta_y);
	  } else if (alpha_min == alpha_y_0) {
	    x = int(q1_x + alpha_min * delta_x);
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
	    if (alpha_x[x + 1] < alpha_y[y + 1] - tol) {
	      alpha_xy[count] = alpha_x[x + 1];
	      x++;
	    } else if (alpha_y[y + 1] < alpha_x[x + 1] - tol) {
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
	  if (alpha_min == alpha_x_0) {
	    x = 0;
	    y = int(q1_y + alpha_min * delta_y);
	  } else if (alpha_min == alpha_y_n) {
	    x = int(q1_x + alpha_min * delta_x);
	    y = ny - 1;
	  } else
	    report_error("something wrong in x+ y-");
	  alpha_xy[count] = alpha_min;
	  i_arr[count] = x;
	  j_arr[count] = y;
	  count++;
	  while (x < nx and y > 0) {
	    i_arr[count] = x;
	    j_arr[count] = y;
	    if (alpha_x[x + 1] < alpha_y[y] - tol) {
	      alpha_xy[count] = alpha_x[x + 1];
	      x++;
	    } else if (alpha_y[y] < alpha_x[x + 1] - tol) {
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
	    y = int(q1_y + alpha_min * delta_y);
	  } else if (alpha_min == alpha_y_0) {
	    x = int(q1_x + alpha_min * delta_x);
	    y = 0;
	  } else
	    report_error("something wrong in x- y+");
	  alpha_xy[count] = alpha_min;
	  i_arr[count] = x;
	  j_arr[count] = y;
	  count++;
	  while (x > 0 and y < ny) {
	    i_arr[count] = x;
	    j_arr[count] = y;
	    if (alpha_x[x] < alpha_y[y + 1] - tol) {
	      alpha_xy[count] = alpha_x[x];
	      x--;
	    } else if (alpha_y[y + 1] < alpha_x[x] - tol) {
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
	    y = int(q1_y + alpha_min * delta_y);
	  } else if (alpha_min == alpha_y_n) {
	    x = int(q1_x + alpha_min * delta_x);
	    y = ny - 1;
	  } else
	    report_error("something wrong in x- y-");
	  alpha_xy[count] = alpha_min;
	  i_arr[count] = x;
	  j_arr[count] = y;
	  count++;
	  while (x > 0 and y > 0) {
	    i_arr[count] = x;
	    j_arr[count] = y;
	    if (alpha_x[x] < alpha_y[y] - tol) {
	      alpha_xy[count] = alpha_x[x];
	      x--;
	    } else if (alpha_y[y] < alpha_x[x] - tol) {
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
#ifdef TESTBP
  if (count > 2 * max_n + 1)
    report_error("forward project overflow");
#endif // TESTBP
  if (count > 0)
    calc_xy_z(pixels, voxels, alpha_xy, i_arr, j_arr, count, a, h, source_z,
	      b_z, d_z, nv, nz, midp, d_conv, delta_z, inv_delz, vox_z);
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
  const real tol = 1e-8;
  if (v_pixels[mid - 1] > source_z - tol and
      v_pixels[mid - 1] < source_z + tol)
    report_error("Oops");
  if (v_pixels[mid] > source_z - tol and
      v_pixels[mid] < source_z + tol)
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

  //#pragma omp parallel for shared(h_pixels, v_pixels, pixels, voxels, angles, d_conv, delta_z, inv_delz, vox_z, voxel_size, grid_offset) firstprivate(n_angles, n_h, n_v, nx_voxels, ny_voxels, nz_voxels) schedule(dynamic)
  for (int a = 0; a < n_angles; a++) {
    /* rotate source and detector positions by current angle */
    real cosa = std::cos(angles[a]);
    real sina = std::sin(angles[a]);
    real p1_x = cosa * source_x - sina * source_y;
    real p1_y = sina * source_x + cosa * source_y;

    for (int h = 0; h < n_h; h++) {
      real p2_x = cosa * detector_x - sina * h_pixels[a];
      real p2_y = sina * detector_x + cosa * h_pixels[a];
      fproject_xy(p1_x, p1_y, p2_x, p2_y, pixels, voxels, grid_offset[0],
		  grid_offset[1], grid_offset[2], voxel_size[0], voxel_size[1],
		  voxel_size[2], nx_voxels, ny_voxels, nz_voxels, a, h,
		  source_z, n_v, mid, d_conv, delta_z, inv_delz, vox_z);
    }
  }
}

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
				const recon_1d &inv_delz, const recon_1d &vox_z)
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
      int k = int(pzbz + alpha_inv * delta_z[v]);
      if (k > 0) {
	recon_type alpha_z = vox_z[k] * inv_delz[v];
	recon_type min_z = std::min(alpha_z, alpha_xy_1[m]);
	voxels[i][j][k] += pixels[a[m]][h[m]][v] * (min_z - alpha_xy_0[m])
	  * d_conv[h[m]][v];
	voxels[i][j][k - 1] += pixels[a[m]][h[m]][v] * (alpha_xy_1[m] - min_z)
	  * d_conv[h[m]][v];
      } else if (k == 0) {
	recon_type alpha_z = vox_z[k] * inv_delz[v];
	recon_type min_z = std::min(alpha_z, alpha_xy_1[m]);
	voxels[i][j][k] += pixels[a[m]][h[m]][v] * (min_z - alpha_xy_0[m])
	  * d_conv[h[m]][v];
      } else
	break;
    }
    for (int v = midp; v < nv; v++) {
      int k = int(pzbz + alpha_inv * delta_z[v]);
      if (k < nzm1) {
	recon_type alpha_z = vox_z[k + 1] * inv_delz[v];
	recon_type min_z = std::min(alpha_z, alpha_xy_1[m]);
	voxels[i][j][k] += pixels[a[m]][h[m]][v] * (min_z - alpha_xy_0[m])
	  * d_conv[h[m]][v];
	voxels[i][j][k + 1] += pixels[a[m]][h[m]][v] * (alpha_xy_1[m] - min_z)
	  * d_conv[h[m]][v];
      } else if (k == nzm1) {
	recon_type alpha_z = vox_z[k + 1] * inv_delz[v];
	recon_type min_z = std::min(alpha_z, alpha_xy_1[m]);
	voxels[i][j][k] += pixels[a[m]][h[m]][v] * (min_z - alpha_xy_0[m])
	  * d_conv[h[m]][v];
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
  so the line representing the detector pixels which in perpendicular to this is
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
  So we solve for t, generate u and apply u to a source line at theta == 0
  in order to map it into a value in the pixel array
  y = - t cphi
*/

/*
  This tries to use the voxel centre and was untested
void bproject_ah(const real source_x, const real source_y,
		 const real detector_x, pixel_data &pixels, voxel_data &voxels,
		 const real b_x, const real b_y, const real b_z,
		 const real d_x, const real d_y, const real d_z,
		 const int nx, const int ny, const int nz, const int i,
		 const int j, const real source_z, const int n_angles,
		 const int n_h, const int n_v, const real_1d &h_pixels,
		 const real_1d &v_pixels, const int midp, const real_1d &cangle,
		 const real_1d &sangle)
{
  const int tol = 1e-8;
  // Todo - check how to calculate this
  const int pix_per_vox = n_v / (nx - 1);
  // How big should the array be
  int count = 0;
  std::vector<int> a_arr(2 * pix_per_vox * n_angles);
  std::vector<int> h_arr(2 * pix_per_vox * n_angles);
  recon_1d alpha_xy_0(2 * pix_per_vox * n_angles);
  recon_1d alpha_xy_1(2 * pix_per_vox * n_angles);
  const real x_0 = b_x + real(i) * d_x;
  const real y_0 = b_y + real(j) * d_y;
  const real x_n = b_x + real(i + 1) * d_x;
  const real y_n = b_y + real(j + 1) * d_y;
  // voxel centre
  real vx = b_x + (real(i) + 0.5) * d_x;
  real vy = b_y + (real(j) + 0.5) * d_y;
  real pixel_step = h_pixels[1] - h_pixels[0];
  for (int a = 0; a < n_angles; a++) {
    real cphi = cangle[a];
    real sphi = sangle[a];
    real px = source_x * cphi - source_y * sphi;
    real py = source_x * sphi + source_y * cphi;
    real qx = detector_x * cphi - source_y * sphi;
    real qy = detector_x * sphi + source_y * cphi;
    real delta_x = vx - px;
    real delta_y = vy - py;
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
    real y = -t * cphi;
    int hmid = int((y - h_pixels[0]) / pixel_step);
    // If it intercepts voxel then calc alpha_xy_0, alpha_xy_1 and store
    // find where line from source through centre of voxel would hit
    // detector and step either side until we miss the voxel
    // Cross-section of voxel at the ray angle versus projected size of steps
    // between pixels should give some idea of size of the side-steps
    // Or can we do something similar with one maximal edge and then
    // just step in one direction? Todo
    // can get some sort of maximal edge by calculating angles of lines
    // through 4 corners and taking lowest (-pi/2 to pi/2 ...?)
    const real p1_x = px;
    const real p1_y = py;
    int pcount = count;
    for (int h = hmid - 1; h >= 0; h--) {
      const real p2_x = cphi * detector_x - sphi * h_pixels[h];
      const real p2_y = sphi * detector_x + cphi * h_pixels[h];
      const real delta_x = p2_x - p1_x;
      const real delta_y = p2_y - p1_y;
      if (std::abs(delta_x) < tol) {
	report_error("Ooops - delta x");
      } else if (std::abs(delta_y) < tol) {
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
	if (alpha_min < alpha_max) {
	  alpha_xy_0[count] = alpha_min;
	  alpha_xy_1[count] = alpha_max;
	  a_arr[count] = a;
	  h_arr[count] = h;
	  count++;
	} else
	  break;
      }
    }
    if (pcount < count) {
      // reorder inserted values so h always increases
      int size = (count - pcount) / 2;
      for (int l = 0; l < size; l++) {
	real tmp = alpha_xy_0[pcount + l];
	alpha_xy_0[pcount + l] = alpha_xy_0[count - l - 1];
	alpha_xy_0[count - l - 1] = tmp;
	tmp = alpha_xy_1[pcount + l];
	alpha_xy_1[pcount + l] = alpha_xy_1[count - l - 1];
	alpha_xy_1[count - l - 1] = tmp;
	int t = h_arr[pcount + l];
	h_arr[pcount + l] = h_arr[count - l - 1];
	h_arr[count - l - 1] = t;
	// a is the same for all
      }
    }
    for (int h = hmid; h < n_h; h++) {
      const real p2_x = cphi * detector_x - sphi * h_pixels[h];
      const real p2_y = sphi * detector_x + cphi * h_pixels[h];
      const real delta_x = p2_x - p1_x;
      const real delta_y = p2_y - p1_y;
      if (std::abs(delta_x) < tol) {
	report_error("Ooops - delta x");
      } else if (std::abs(delta_y) < tol) {
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
	if (alpha_min < alpha_max) {
	  alpha_xy_0[count] = alpha_min;
	  alpha_xy_1[count] = alpha_max;
	  a_arr[count] = a;
	  h_arr[count] = h;
	  count++;
	} else
	  break;
      }
    }
  }
#ifdef TESTBP
  if (count > 2 * pix_per_vox * n_angles)
    report_error("back project overflow");
#endif // TESTBP
  if (count > 0)
    do_z();
}
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
  const int tol = 1e-8;
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
    real y00 = -t * cphi;
    int h00 = int((y00 - h_pixels[0]) / pixel_step);
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
    real y01 = -t * cphi;
    int h01 = int((y01 - h_pixels[0]) / pixel_step);
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
    real y10 = -t * cphi;
    int h10 = int((y10 - h_pixels[0]) / pixel_step);
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
    real y11 = -t * cphi;
    int h11 = int((y11 - h_pixels[0]) / pixel_step);
    hmin = std::max(std::min(hmin, h11), 0);
    hmax = std::min(std::max(hmax, h11), n_h - 1);
    //int hmin = std::max(std::min(std::min(h00, h01), std::min(h10, h11)), 0);
    //int hmax = std::min(std::max(std::max(h00, h01), std::max(h10, h11)),
    //		n_h - 1);
    // If it intercepts voxel then calc alpha_xy_0, alpha_xy_1 and store
    const real p1_x = px;
    const real p1_y = py;
    for (int h = hmin; h <= hmax; h++) {
      const real p2_x = cphi * detector_x - sphi * h_pixels[h];
      const real p2_y = sphi * detector_x + cphi * h_pixels[h];
      const real delta_x = p2_x - p1_x;
      const real delta_y = p2_y - p1_y;
      if (std::abs(delta_x) < tol) {
	report_error("Ooops - delta x");
      } else if (std::abs(delta_y) < tol) {
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
	if (alpha_min < alpha_max) {
	  alpha_xy_0[count] = alpha_min;
	  alpha_xy_1[count] = alpha_max;
	  a_arr[count] = a;
	  h_arr[count] = h;
	  count++;
	} else
	  break;
      }
    }
  }
#ifdef TESTBP
  if (count > 2 * pix_per_vox * n_angles)
    report_error("back project overflow");
#endif // TESTBP
  if (count > 0)
    calc_ah_z(pixels, voxels, alpha_xy_0, alpha_xy_1, a_arr, h_arr, count,
	      i, j, source_z, b_z, d_z, n_v, nz, midp, d_conv, delta_z,
	      inv_delz, vox_z);
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
