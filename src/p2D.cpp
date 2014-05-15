
// Ugly copy of projection algorithm, could simplify to 2D
// Todo - how to avoid a copy?
#include <map>
#include <vector>
#include <cmath>
#ifdef MATLAB_MEX_FILE
#  include "mex_types.hpp"
#else
#  include "base_types.hpp"
#endif // mex
#include "instruments.hpp"
#include "ui_calls.hpp"

extern void test_voxel(const int i, const int j, const int k,
		       const int a, const int v, const int h, const double d,
		       const double l, const double voxel_origin[3],
		       const double voxel_size[3], const double v_pos,
		       const double h_pos, const double cphi,
		       const double sphi, const int nx, const int ny,
		       const int nz, const char message[]);

void CCPi::parallel_beam::f2D(const real_1d &h_pixels, const real_1d &v_pixels,
			      const real_1d &angles, const int n_angles,
			      const int nh_pixels, const int nv_pixels,
			      const real vox_origin[3], const real vox_size[3],
			      const int nx, const int ny, const int nz,
			      pixel_data &pixels, voxel_data &voxels)
{
  //sl_int nyz = sl_int(ny) * sl_int(nz);
  const int pixels_per_voxel = nv_pixels / nz;
  std::vector<real> cangles(n_angles);
  for (int i = 0; i < n_angles; i++)
    cangles[i] = std::cos(angles[i]);
  std::vector<real> sangles(n_angles);
  for (int i = 0; i < n_angles; i++)
    sangles[i] = std::sin(angles[i]);
  //real inv_pixel_step = 1.0 / (h_pixels[1] - h_pixels[0]);
  std::vector<real> vox_lengths(n_angles);
  for (int i = 0; i < n_angles; i++) {
    if (std::abs(cangles[i]) > std::abs(sangles[i]))
      vox_lengths[i] = vox_size[0] / std::abs(cangles[i]);
    else
      vox_lengths[i] = vox_size[0] / std::abs(sangles[i]);
  }
  real x_0 = vox_origin[0];
  real x_n = vox_origin[0] + real(nx) * vox_size[0];
  real y_0 = vox_origin[1];
  real y_n = vox_origin[1] + real(ny) * vox_size[1];
  // assume square voxels
  real inv_vox_step = 1.0 / vox_size[0];
  // based on Xaolin Wu's antialised Bresenhams algorithm
#pragma omp parallel for shared(h_pixels, pixels, cangles, sangles, vox_lengths) firstprivate(pixels_per_voxel, nx, ny, nz, nh_pixels, inv_vox_step, x_0, x_n, y_0, y_n) schedule(dynamic)
  for (int a = 0; a < n_angles; a++) {
#ifdef TESTBP
    std::vector<real> lengths(2 * nx);
    std::vector<int> idex(2 * nx);
    std::vector<int> jdex(2 * nx);
#endif // TESTBP
    const real tol = 1e-11; // needs to be > nx * epsilon
    real cphi = cangles[a];
    real sphi = sangles[a];
    real L = vox_lengths[a];
    bool steep = std::abs(sphi) > std::abs(cphi);
    for (int h = 0; h < nh_pixels; h++) {
      //sl_int pix_av = (sl_int(a) * sl_int(nh_pixels) + sl_int(h)) * sl_int(nv_pixels);
      //find end points of line (x0,y0), (x1,y1), scaled from 0 to n, step 1.0
      // x = -sphi * hp[h] + r cphi, y = cphi * hp[h] + r * sphi
      real x0;
      real x1;
      real y0;
      real y1;
      if (sphi < 0.0) { // 180-360
	if (cphi < 0.0) { // 180-270
	  if (steep) {
	    y0 = y_n - tol;
	    real r0 = (y_n - cphi * h_pixels[h]) / sphi;
	    x0 = -sphi * h_pixels[h] + r0 * cphi;
	    if (x0 >= x_n and std::abs(cphi) > 1e-8) {
	      x0 = x_n - tol;
	      r0 = (x_n + sphi * h_pixels[h]) / cphi;
	      y0 = cphi * h_pixels[h] + r0 * sphi;
	    }
	    y1 = y_0;
	    real r1 = (y_0 - cphi * h_pixels[h]) / sphi;
	    x1 = -sphi * h_pixels[h] + r1 * cphi;
	    if (x1 < x_0 and std::abs(cphi) > 1e-8) {
	      x1 = x_0;
	      r1 = (x_0 + sphi * h_pixels[h]) / cphi;
	      y1 = cphi * h_pixels[h] + r1 * sphi;
	    }
	  } else {
	    x0 = x_n - tol;
	    real r0 = (x_n + sphi * h_pixels[h]) / cphi;
	    y0 = cphi * h_pixels[h] + r0 * sphi;
	    if (y0 >= y_n and std::abs(sphi) > 1e-8) {
	      y0 = y_n - tol;
	      r0 = (y_n - cphi * h_pixels[h]) / sphi;
	      x0 = -sphi * h_pixels[h] + r0 * cphi;
	    }
	    x1 = x_0;
	    real r1 = (x_0 + sphi * h_pixels[h]) / cphi;
	    y1 = cphi * h_pixels[h] + r1 * sphi;
	    if (y1 < y_0 and std::abs(sphi) > 1e-8) {
	      y1 = y_0;
	      r1 = (y_0 - cphi * h_pixels[h]) / sphi;
	      x1 = -sphi * h_pixels[h] + r1 * cphi;
	    }
	  }
	} else { // 270-360
	  if (steep) {
	    y0 = y_n - tol;
	    real r0 = (y_n - cphi * h_pixels[h]) / sphi;
	    x0 = -sphi * h_pixels[h] + r0 * cphi;
	    if (x0 < x_0 and std::abs(cphi) > 1e-8) {
	      x0 = x_0;
	      r0 = (x_0 + sphi * h_pixels[h]) / cphi;
	      y0 = cphi * h_pixels[h] + r0 * sphi;
	    }
	    y1 = y_0;
	    real r1 = (y_0 - cphi * h_pixels[h]) / sphi;
	    x1 = -sphi * h_pixels[h] + r1 * cphi;
	    if (x1 >= x_n and std::abs(cphi) > 1e-8) {
	      x1 = x_n - tol;
	      r1 = (x_n + sphi * h_pixels[h]) / cphi;
	      y1 = cphi * h_pixels[h] + r1 * sphi;
	    }
	  } else {
	    x0 = x_0;
	    real r0 = (x_0 + sphi * h_pixels[h]) / cphi;
	    y0 = cphi * h_pixels[h] + r0 * sphi;
	    if (y0 >= y_n and std::abs(sphi) > 1e-8) {
	      y0 = y_n - tol;
	      r0 = (y_n - cphi * h_pixels[h]) / sphi;
	      x0 = -sphi * h_pixels[h] + r0 * cphi;
	    }
	    x1 = x_n - tol;
	    real r1 = (x_n + sphi * h_pixels[h]) / cphi;
	    y1 = cphi * h_pixels[h] + r1 * sphi;
	    if (y1 < y_0 and std::abs(sphi) > 1e-8) {
	      y1 = y_0;
	      r1 = (y_0 - cphi * h_pixels[h]) / sphi;
	      x1 = -sphi * h_pixels[h] + r1 * cphi;
	    }
	  }
	}
      } else {
	if (cphi < 0.0) { // 90-180
	  if (steep) {
	    y0 = y_0;
	    real r0 = (y_0 - cphi * h_pixels[h]) / sphi;
	    x0 = -sphi * h_pixels[h] + r0 * cphi;
	    if (x0 >= x_n and std::abs(cphi) > 1e-8) {
	      x0 = x_n - tol;
	      r0 = (x_n + sphi * h_pixels[h]) / cphi;
	      y0 = cphi * h_pixels[h] + r0 * sphi;
	    }
	    y1 = y_n - tol;
	    real r1 = (y_n - cphi * h_pixels[h]) / sphi;
	    x1 = -sphi * h_pixels[h] + r1 * cphi;
	    if (x1 < x_0 and std::abs(cphi) > 1e-8) {
	      x1 = x_0;
	      r1 = (x_0 + sphi * h_pixels[h]) / cphi;
	      y1 = cphi * h_pixels[h] + r1 * sphi;
	    }
	  } else {
	    x0 = x_n - tol;
	    real r0 = (x_n + sphi * h_pixels[h]) / cphi;
	    y0 = cphi * h_pixels[h] + r0 * sphi;
	    if (y0 < y_0 and std::abs(sphi) > 1e-8) {
	      y0 = y_0;
	      r0 = (y_0 - cphi * h_pixels[h]) / sphi;
	      x0 = -sphi * h_pixels[h] + r0 * cphi;
	    }
	    x1 = x_0;
	    real r1 = (x_0 + sphi * h_pixels[h]) / cphi;
	    y1 = cphi * h_pixels[h] + r1 * sphi;
	    if (y1 >= y_n and std::abs(sphi) > 1e-8) {
	      y1 = y_n - tol;
	      r1 = (y_n - cphi * h_pixels[h]) / sphi;
	      x1 = -sphi * h_pixels[h] + r1 * cphi;
	    }
	  }
	} else { // 0-90
	  if (steep) {
	    y0 = y_0;
	    real r0 = (y_0 - cphi * h_pixels[h]) / sphi;
	    x0 = -sphi * h_pixels[h] + r0 * cphi;
	    if (x0 < x_0 and std::abs(cphi) > 1e-8) {
	      x0 = x_0;
	      r0 = (x_0 + sphi * h_pixels[h]) / cphi;
	      y0 = cphi * h_pixels[h] + r0 * sphi;
	    }
	    y1 = y_n - tol;
	    real r1 = (y_n - cphi * h_pixels[h]) / sphi;
	    x1 = -sphi * h_pixels[h] + r1 * cphi;
	    if (x1 >= x_n and std::abs(cphi) > 1e-8) {
	      x1 = x_n - tol;
	      r1 = (x_n + sphi * h_pixels[h]) / cphi;
	      y1 = cphi * h_pixels[h] + r1 * sphi;
	    }
	  } else {
	    x0 = x_0;
	    real r0 = (x_0 + sphi * h_pixels[h]) / cphi;
	    y0 = cphi * h_pixels[h] + r0 * sphi;
	    if (y0 < y_0 and std::abs(sphi) > 1e-8) {
	      y0 = y_0;
	      r0 = (y_0 - cphi * h_pixels[h]) / sphi;
	      x0 = -sphi * h_pixels[h] + r0 * cphi;
	    }
	    x1 = x_n - tol;
	    real r1 = (x_n + sphi * h_pixels[h]) / cphi;
	    y1 = cphi * h_pixels[h] + r1 * sphi;
	    if (y1 >= y_n and std::abs(sphi) > 1e-8) {
	      y1 = y_n - tol;
	      r1 = (y_n - cphi * h_pixels[h]) / sphi;
	      x1 = -sphi * h_pixels[h] + r1 * cphi;
	    }
	  }
	}
      }
      //fprintf(stderr, "%f %f %f %f, %f %f %f %f\n", x0, y0, x1, y1,
      //      x_0, y_0, x_n, y_n);
      if (x0 >= x_0 and x1 < x_n and y0 >= y_0 and y1 < y_n) {
	// rescale from 0 to n, step 1, shifted to 0,0
	x0 = (x0 - x_0) * inv_vox_step;
	x1 = (x1 - x_0) * inv_vox_step;
	y0 = (y0 - y_0) * inv_vox_step;
	y1 = (y1 - y_0) * inv_vox_step;
	// Now start the algorithm
	//bool steep = std::abs(y1 - y0) > std::abs(x1 - x0);
	if (steep) {
	  real t = x0;
	  x0 = y0;
	  y0  = t;
	  t = x1;
	  x1 = y1;
	  y1 = t;
	}
	if (x0 > x1) {
	  real t = x0;
	  x0 = x1;
	  x1 = t;
	  t = y0;
	  y0 = y1;
	  y1 = t;
	}
	real dx = x1 - x0;
	real dy = y1 - y0;
#ifdef TESTBP
	int cnt = 0;
#endif // TESTBP
	real gradient = dy / dx;
	// first endpoint
	real xend = std::floor(x0);
	//real yend = y0 + gradient * (xend - x0);
	int xpxl1 = int(x0);
	int y = int(y0);
	real xgap = x0 - xend;
	real intery = y0 - std::floor(y0);
	const real gtol = 5e-9;
	if (xgap > gtol) {
	  xgap = 1.0 - xgap;
	  real len = xgap * L;
	  if (steep) {
	    //use(y, xpxl1, L * xgap);
	    sl_int k_range = 0; //pix_av;
	    //sl_int vox_y = sl_int(y) * nyz + sl_int(xpxl1) * sl_int(nz);
	    for (int k = 0; k < nz; k++) {
	      for (int v = 0; v < pixels_per_voxel; v++)
		pixels[a][h][k_range + v] += voxels[y][xpxl1][k] * len;
	      k_range += pixels_per_voxel;
	    }
#ifdef TESTBP
	    lengths[cnt] = len;
	    idex[cnt] = y;
	    jdex[cnt] = xpxl1;
	    cnt++;
#endif // TESTBP
	  } else {
	    //use(xpxl1, y, L * xgap);
	    sl_int k_range = 0; //pix_av;
	    //sl_int vox_y = sl_int(xpxl1) * nyz + sl_int(y) * sl_int(nz);
	    for (int k = 0; k < nz; k++) {
	      for (int v = 0; v < pixels_per_voxel; v++)
		pixels[a][h][k_range + v] += voxels[xpxl1][y][k] * len;
	      k_range += pixels_per_voxel;
	    }
#ifdef TESTBP
	    lengths[cnt] = len;
	    idex[cnt] = xpxl1;
	    jdex[cnt] = y;
	    cnt++;
#endif // TESTBP
	  }
	  intery += gradient * xgap;
	  if (gradient < 0.0) {
	    if (intery >= 1.0) {
	      y++;
	      intery -= 1.0;
	    }
	  } else {
	    if (intery < 0.0) {
	      y--;
	      intery += 1.0;
	    }
	  }
	  xpxl1++;
	}
	int xpxl2 = int(x1);
	xend = std::floor(x1);
	//yend = y + gradient * (x1 - xend);
	xgap = x1 - xend;
	bool extra = true;
	if (xgap > 1.0 - gtol) {
	  xpxl2++;
	  extra = false;
	}
	//int ypxl2 = int(y1);
	// main loop
	if (gradient < 0.0) {
	  real Lscale = L / -gradient;
	  for (int x = xpxl1; x < xpxl2; x++) {
	    real under = intery;
	    intery += gradient;
	    if (intery < 0.0) {
	      real over = -gradient - under;
	      intery += 1.0;
	      if (steep) {
		//use(y, x, Lscale * under);
		real len = Lscale * under;
		sl_int k_range = 0; //pix_av;
		//sl_int vox_y = sl_int(y) * nyz + sl_int(x) * sl_int(nz);
		for (int k = 0; k < nz; k++) {
		  for (int v = 0; v < pixels_per_voxel; v++)
		    pixels[a][h][k_range + v] += voxels[y][x][k] * len;
		  k_range += pixels_per_voxel;
		}
#ifdef TESTBP
		lengths[cnt] = len;
		idex[cnt] = y;
		jdex[cnt] = x;
		cnt++;
#endif // TESTBP
		// L - L * under / gradient;
		//use(y + 1, x, Lscale * over);
		len = Lscale * over;
		k_range = 0; //pix_av;
		//vox_y -= nyz;
		for (int k = 0; k < nz; k++) {
		  for (int v = 0; v < pixels_per_voxel; v++)
		    pixels[a][h][k_range + v] += voxels[y + 1][x][k] * len;
		  k_range += pixels_per_voxel;
		}
#ifdef TESTBP
		lengths[cnt] = len;
		idex[cnt] = y - 1;
		jdex[cnt] = x;
		cnt++;
#endif // TESTBP
	      } else {
		//use(x, y, Lscale * under);
		real len = Lscale * under;
		sl_int k_range = 0; //pix_av;
		//sl_int vox_y = sl_int(x) * nyz + sl_int(y) * sl_int(nz);
		for (int k = 0; k < nz; k++) {
		  for (int v = 0; v < pixels_per_voxel; v++)
		    pixels[a][h][k_range + v] += voxels[x][y][k] * len;
		  k_range += pixels_per_voxel;
		}
#ifdef TESTBP
		lengths[cnt] = len;
		idex[cnt] = x;
		jdex[cnt] = y;
		cnt++;
#endif // TESTBP
		// L - L * under / gradient;
		//use(x, y + 1, Lscale * over);
		len = Lscale * over;
		k_range = 0; //pix_av;
		//vox_y -= sl_int(nz);
		for (int k = 0; k < nz; k++) {
		  for (int v = 0; v < pixels_per_voxel; v++)
		    pixels[a][h][k_range + v] += voxels[x][y + 1][k] * len;
		  k_range += pixels_per_voxel;
		}
#ifdef TESTBP
		lengths[cnt] = len;
		idex[cnt] = x;
		jdex[cnt] = y - 1;
		cnt++;
#endif // TESTBP
	      }
	      y--;
	    } else {
	      real len = L;
	      if (steep) {
		//use(y, x, L);
		sl_int k_range = 0; //pix_av;
		//sl_int vox_y = sl_int(y) * nyz + sl_int(x) * sl_int(nz);
		for (int k = 0; k < nz; k++) {
		  for (int v = 0; v < pixels_per_voxel; v++)
		    pixels[a][h][k_range + v] += voxels[y][x][k] * len;
		  k_range += pixels_per_voxel;
		}
#ifdef TESTBP
		lengths[cnt] = len;
		idex[cnt] = y;
		jdex[cnt] = x;
		cnt++;
#endif // TESTBP
	      } else {
		//use(x, y, L);
		sl_int k_range = 0; //pix_av;
		//sl_int vox_y = sl_int(x) * nyz + sl_int(y) * sl_int(nz);
		for (int k = 0; k < nz; k++) {
		  for (int v = 0; v < pixels_per_voxel; v++)
		    pixels[a][h][k_range + v] += voxels[x][y][k] * len;
		  k_range += pixels_per_voxel;
		}
#ifdef TESTBP
		lengths[cnt] = len;
		idex[cnt] = x;
		jdex[cnt] = y;
		cnt++;
#endif // TESTBP
	      }
	    }
	  }
	} else {
	  real Lscale = L / gradient;
	  for (int x = xpxl1; x < xpxl2; x++) {
	    real under = 1.0 - intery;
	    intery += gradient;
	    if (intery >= 1.0) {
	      real over = gradient - under;
	      intery -= 1.0;
	      if (steep) {
		//use(y, x, Lscale * under);
		real len = Lscale * under;
		sl_int k_range = 0; //pix_av;
		//sl_int vox_y = sl_int(y) * nyz + sl_int(x) * sl_int(nz);
		for (int k = 0; k < nz; k++) {
		  for (int v = 0; v < pixels_per_voxel; v++)
		    pixels[a][h][k_range + v] += voxels[y][x][k] * len;
		  k_range += pixels_per_voxel;
		}
#ifdef TESTBP
		lengths[cnt] = len;
		idex[cnt] = y;
		jdex[cnt] = x;
		cnt++;
#endif // TESTBP
		// L - L * under / gradient;
		//use(y + 1, x, Lscale * over);
		len = Lscale * over;
		k_range = 0; //pix_av;
		//vox_y += nyz;
		for (int k = 0; k < nz; k++) {
		  for (int v = 0; v < pixels_per_voxel; v++)
		    pixels[a][h][k_range + v] += voxels[y + 1][x][k] * len;
		  k_range += pixels_per_voxel;
		}
#ifdef TESTBP
		lengths[cnt] = len;
		idex[cnt] = y + 1;
		jdex[cnt] = x;
		cnt++;
#endif // TESTBP
	      } else {
		//use(x, y, Lscale * under);
		real len = Lscale * under;
		sl_int k_range = 0; //pix_av;
		//sl_int vox_y = sl_int(x) * nyz + sl_int(y) * sl_int(nz);
		for (int k = 0; k < nz; k++) {
		  for (int v = 0; v < pixels_per_voxel; v++)
		    pixels[a][h][k_range + v] += voxels[x][y][k] * len;
		  k_range += pixels_per_voxel;
		}
#ifdef TESTBP
		lengths[cnt] = len;
		idex[cnt] = x;
		jdex[cnt] = y;
		cnt++;
#endif // TESTBP
		// L - L * under / gradient;
		//use(x, y + 1, Lscale * over);
		len = Lscale * over;
		k_range = 0; //pix_av;
		//vox_y += sl_int(nz);
		for (int k = 0; k < nz; k++) {
		  for (int v = 0; v < pixels_per_voxel; v++)
		    pixels[a][h][k_range + v] += voxels[x][y + 1][k] * len;
		  k_range += pixels_per_voxel;
		}
#ifdef TESTBP
		lengths[cnt] = len;
		idex[cnt] = x;
		jdex[cnt] = y + 1;
		cnt++;
#endif // TESTBP
	      }
	      y++;
	    } else {
	      real len = L;
	      if (steep) {
		//use(y, x, L);
		sl_int k_range = 0; //pix_av;
		//sl_int vox_y = sl_int(y) * nyz + sl_int(x) * sl_int(nz);
		for (int k = 0; k < nz; k++) {
		  for (int v = 0; v < pixels_per_voxel; v++)
		    pixels[a][h][k_range + v] += voxels[y][x][k] * len;
		  k_range += pixels_per_voxel;
		}
#ifdef TESTBP
		lengths[cnt] = len;
		idex[cnt] = y;
		jdex[cnt] = x;
		cnt++;
#endif // TESTBP
	      } else {
		//use(x, y, L);
		sl_int k_range = 0; //pix_av;
		//sl_int vox_y = sl_int(x) * nyz + sl_int(y) * sl_int(nz);
		for (int k = 0; k < nz; k++) {
		  for (int v = 0; v < pixels_per_voxel; v++)
		    pixels[a][h][k_range + v] += voxels[x][y][k] * len;
		  k_range += pixels_per_voxel;
		}
#ifdef TESTBP
		lengths[cnt] = len;
		idex[cnt] = x;
		jdex[cnt] = y;
		cnt++;
#endif // TESTBP
	      }
	    }
	  }
	}
	if (extra and xpxl2 > xpxl1) {
	  xend = std::floor(x1);
	  //yend = y + gradient * (x1 - xend);
	  xgap = x1 - xend;
	  // partial x step to y edge
	  real len = L * xgap;
	  if (steep) {
	    //use(y, xpxl2, L * xgap);
	    sl_int k_range = 0; //pix_av;
	    //sl_int vox_y = sl_int(y) * nyz + sl_int(xpxl2) * sl_int(nz);
	    for (int k = 0; k < nz; k++) {
	      for (int v = 0; v < pixels_per_voxel; v++)
		pixels[a][h][k_range + v] += voxels[y][xpxl2][k] * len;
	      k_range += pixels_per_voxel;
	    }
#ifdef TESTBP
	    lengths[cnt] = len;
	    idex[cnt] = y;
	    jdex[cnt] = xpxl2;
	    cnt++;
#endif // TESTBP
	  } else {
	    //use(xpxl2, y, L * xgap);
	    sl_int k_range = 0; //pix_av;
	    //sl_int vox_y = sl_int(xpxl2) * nyz + sl_int(y) * sl_int(nz);
	    for (int k = 0; k < nz; k++) {
	      for (int v = 0; v < pixels_per_voxel; v++)
		pixels[a][h][k_range + v] += voxels[xpxl2][y][k] * len;
	      k_range += pixels_per_voxel;
	    }
#ifdef TESTBP
	    lengths[cnt] = len;
	    idex[cnt] = xpxl2;
	    jdex[cnt] = y;
	    cnt++;
#endif // TESTBP
	  }
	}
#ifdef TESTBP
	for (int c = 0; c < cnt; c++) {
	  test_voxel(idex[c], jdex[c], 0, a, 0, h, 1.0, lengths[c],
		     vox_origin, vox_size, v_pixels[0], h_pixels[h],
		     cphi, sphi, nx, ny, nz, "xy");
	  if (idex[c] - 1 >= 0) {
	    int d;
	    for (d = 0; d < cnt; d++) {
	      if (idex[d] == idex[c] - 1 and jdex[d] == jdex[c])
		break;
	    }
	    if (d == cnt)
	      test_voxel(idex[c]-1, jdex[c], 0, a, 0, h, 1.0, 0.0,
			 vox_origin, vox_size, v_pixels[0], h_pixels[h],
			 cphi, sphi, nx, ny, nz, "i-");
	  }
	  if (idex[c] + 1 < nx) {
	    int d;
	    for (d = 0; d < cnt; d++) {
	      if (idex[d] == idex[c] + 1 and jdex[d] == jdex[c])
		break;
	    }
	    if (d == cnt)
	      test_voxel(idex[c]+1, jdex[c], 0, a, 0, h, 1.0, 0.0,
			 vox_origin, vox_size, v_pixels[0], h_pixels[h],
			 cphi, sphi, nx, ny, nz, "i+");
	  }
	  if (jdex[c] - 1 >= 0) {
	    int d;
	    for (d = 0; d < cnt; d++) {
	      if (idex[d] == idex[c] and jdex[d] == jdex[c] - 1)
		break;
	    }
	    if (d == cnt)
	      test_voxel(idex[c], jdex[c]-1, 0, a, 0, h, 1.0, 0.0,
			 vox_origin, vox_size, v_pixels[0], h_pixels[h],
			 cphi, sphi, nx, ny, nz, "j-");
	  }
	  if (jdex[c] + 1 < ny) {
	    int d;
	    for (d = 0; d < cnt; d++) {
	      if (idex[d] == idex[c] and jdex[d] == jdex[c] + 1)
		break;
	    }
	    if (d == cnt)
	      test_voxel(idex[c], jdex[c]+1, 0, a, 0, h, 1.0, 0.0,
			 vox_origin, vox_size, v_pixels[0], h_pixels[h],
			 cphi, sphi, nx, ny, nz, "j+");
	  }
	}
#endif // TESTBP
      }
    }
  }
}

// Todo - think about scaling to 1.0 ideas, loop blocking for cache reuse?
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
  //const real tol = 1e-13;
  //sl_int nyz = sl_int(ny) * sl_int(nz);
  const int pixels_per_voxel = nv_pixels / nz;
  std::vector<real> cangles(n_angles);
  for (int i = 0; i < n_angles; i++)
    cangles[i] = std::cos(angles[i]);
  std::vector<real> sangles(n_angles);
  for (int i = 0; i < n_angles; i++)
    sangles[i] = std::sin(angles[i]);
  // centred positions
  std::vector<real> x0_coords(nx + 1);
  for (int i = 0; i <= nx; i++)
    x0_coords[i] = vox_origin[0] + real(i) * vox_size[0];
  std::vector<real> y0_coords(ny + 1);
  for (int i = 0; i <= ny; i++)
    y0_coords[i] = vox_origin[1] + real(i) * vox_size[1];
  std::vector<real> z_coords(nz + 1);
  for (int i = 0; i <= nz; i++)
    z_coords[i] = vox_origin[2] + real(i) * vox_size[2];
  real inv_pixel_step = 1.0 / (h_pixels[1] - h_pixels[0]);
  std::vector<real> vox_lengths(n_angles);
  for (int i = 0; i < n_angles; i++) {
    if (std::abs(cangles[i]) > std::abs(sangles[i]))
      vox_lengths[i] = vox_size[0] / std::abs(cangles[i]);
    else
      vox_lengths[i] = vox_size[0] / std::abs(sangles[i]);
  }
  std::vector<real> v_coords(nv_pixels + 1);
  for (int i = 0; i < nv_pixels; i++)
    v_coords[i] = v_pixels[i];
  v_coords[nv_pixels] = z_coords[nz] + 1.0;
  if (v_pixels[0] < z_coords[0] or v_pixels[pixels_per_voxel] < z_coords[1])
    report_error("Oops - bad back projection");
#pragma omp parallel for shared(h_pixels, pixels, x0_coords, y0_coords, cangles, sangles, vox_lengths) firstprivate(pixels_per_voxel, nx, ny, nz, nh_pixels, inv_pixel_step) schedule(dynamic)
  for (int i = 0; i < nx; i++) {
#ifdef TESTBP
    std::vector<real> lengths(nh_pixels);
#endif // TESTBP
    //sl_int vox_xy = i * nyz;
    for (int j = 0; j < ny; j++) {
      //sl_int vox_y = vox_xy + j * sl_int(nz);
      for (int a = 0; a < n_angles; a++) {
	//sl_int pix_av = a * sl_int(nv_pixels) * sl_int(nh_pixels);
	real cphi = cangles[a];
	real sphi = sangles[a];
	real L = vox_lengths[a];
	real xx_0 = - x0_coords[i] * sphi;
	real xx_1 = - x0_coords[i + 1] * sphi;
	real x0_0 = - x0_coords[0] * sphi;
	real yy_0 = y0_coords[j] * cphi;
	real yy_1 = y0_coords[j + 1] * cphi;
	real y_00 = xx_0 + yy_0;
	real y_01 = xx_0 + yy_1;
	real y_10 = xx_1 + yy_0;
	real y_11 = xx_1 + yy_1;
	real y_top = 0.0;
	real y_l1;
	real y_l2;
	real y_bot = 0.0;
	real yj = 0.0;
	if (sphi >= 0.0) {
	  // <= 180, scanning from 0 to nx decreases p
	  // orientation of voxels is the same for all so can get rotated shape
	  // from any single one - use 0.0
	  //int jm;
	  if (cphi >= 0.0) { // 0-90
	    //jm = 1;
	    yj = yy_1;
	    if (cphi >= sphi) { // 0-45
	      y_l1 = y_01 - y_11;
	      y_l2 = y_01 - y_00;
	      y_bot = y_01 - y_10;
	    } else {
	      y_l1 = y_01 - y_00;
	      y_l2 = y_01 - y_11;
	      y_bot = y_01 - y_10;
	    }
	  } else {
	    //jm = 0;
	    yj = yy_0;
	    if (std::abs(cphi) > sphi) { // 135-180
	      y_l1 = y_00 - y_10;
	      y_l2 = y_00 - y_01;
	      y_bot = y_00 - y_11;
	    } else {
	      y_l1 = y_00 - y_01;
	      y_l2 = y_00 - y_10;
	      y_bot = y_00 - y_11;
	    }
	  }
	  real inv_y_l1 = L / y_l1;
	  int p = int((x0_0 + yj - h_pixels[0]) * inv_pixel_step);
	  if (p >= nh_pixels)
	    p = nh_pixels - 1;
	  else if (p < 0)
	    continue;
	  y_top = xx_0 + yj;
	  while (h_pixels[p] > y_top) {
	    p--;
	    if (p < 0)
	      break;
	  }
	  if (p < 0)
	    break;
	  real yb = y_top - y_bot;
	  real y1 = y_top - y_l1;
	  real y2 = y_top - y_l2;
#ifdef TESTBP
	  int cnt = 0;
#endif // TESTBP
	  for (int q = p; q > -1; q--) {
	    if (h_pixels[q] < yb)
	      break;
	    else {
	      real len;
	      if (h_pixels[q] < y2)
		len = (h_pixels[q] - yb) * inv_y_l1;
	      else if (h_pixels[q] < y1)
		len = L;
	      else
		len = (y_top - h_pixels[q]) * inv_y_l1;
#ifdef TESTBP
	      lengths[cnt] = len;
#endif // TESTBP
	      sl_int k_range = 0; //pix_av + q * sl_int(nv_pixels);
	      for (int k = 0; k < nz; k++) {
		for (int v = 0; v < pixels_per_voxel; v++)
		  voxels[i][j][k] += pixels[a][q][k_range + v] * len;
		k_range += pixels_per_voxel;
	      }
#ifdef TESTBP
	      cnt++;
#endif // TESTBP
	    }
	  }
#ifdef TESTBP
	  if (p < nh_pixels - 1)
	    test_voxel(i, j, 0, a, 0, p+1, 0.0, 0.0,
		       vox_origin, vox_size, v_coords[0], h_pixels[p+1],
		       cphi, sphi, nx, ny, nz, "p+");
	  for (int q = 0; q < cnt; q++) {
	    test_voxel(i, j, 0, a, 0, p-q, 1.0, lengths[q],
		       vox_origin, vox_size, v_coords[0], h_pixels[p-q],
		       cphi, sphi, nx, ny, nz, "pq");
	  }
	  if (p-cnt >= 0) {
	    test_voxel(i, j, 0, a, 0, p-cnt, 0.0, 0.0,
		       vox_origin, vox_size, v_coords[0], h_pixels[p-cnt],
		       cphi, sphi, nx, ny, nz, "p-");
	  }
#endif // TESTBP
	} else {
	  // 0-nx increases p
	  //int jm;
	  if (cphi < 0.0) { // 180-270
	    //jm = 1;
	    yj = yy_1;
	    if (cphi < sphi) { // 180-225
	      y_l1 = y_11 - y_01;
	      y_l2 = y_00 - y_01;
	      y_top = y_10 - y_01;
	    } else {
	      y_l1 = y_00 - y_01;
	      y_l2 = y_11 - y_01;
	      y_top = y_10 - y_01;
	    }
	  } else { //270-360
	    //jm = 0;
	    yj = yy_0;
	    if (cphi > std::abs(sphi)) { // 315-360
	      y_l1 = y_10 - y_00;
	      y_l2 = y_01 - y_00;
	      y_top = y_11 - y_00;
	    } else {
	      y_l1 = y_01 - y_00;
	      y_l2 = y_10 - y_00;
	      y_top = y_11 - y_00;
	    }
	  }
	  real inv_y_l1 = L / y_l1;
	  int p = int((x0_0 + yj - h_pixels[0]) * inv_pixel_step);
	  if (p < 0)
	    p = 0;
	  else if (p >= nh_pixels)
	    continue;
	  y_bot = xx_0 + yj;
	  while (h_pixels[p] < y_bot) {
	    p++;
	    if (p >= nh_pixels)
	      break;
	  }
	  if (p >= nh_pixels)
	    break;
	  real yt = y_bot + y_top;
	  real y1 = y_bot + y_l1;
	  real y2 = y_bot + y_l2;
#ifdef TESTBP
	  int cnt = 0;
#endif // TESTBP
	  for (int q = p; q < nh_pixels; q++) {
	    if (h_pixels[q] >= yt)
	      break;
	    else {
	      real len;
	      if (h_pixels[q] < y1)
		len = (h_pixels[q] - y_bot) * inv_y_l1;
	      else if (h_pixels[q] < y2)
		len = L;
	      else
		len = (yt - h_pixels[q]) * inv_y_l1;
#ifdef TESTBP
	      lengths[cnt] = len;
#endif // TESTBP
	      sl_int k_range = 0; //pix_av + q * sl_int(nv_pixels);
	      for (int k = 0; k < nz; k++) {
		for (int v = 0; v < pixels_per_voxel; v++)
		  voxels[i][j][k] += pixels[a][q][k_range + v] * len;
		k_range += pixels_per_voxel;
	      }
#ifdef TESTBP
	      cnt++;
#endif // TESTBP
	    }
	  }
#ifdef TESTBP
	  if (p > 0)
	    test_voxel(i, j, 0, a, 0, p-1, 0.0, 0.0,
		       vox_origin, vox_size, v_coords[0], h_pixels[p-1],
		       cphi, sphi, nx, ny, nz, "p-");
	  for (int q = 0; q < cnt; q++) {
	    test_voxel(i, j, 0, a, 0, p+q, 1.0, lengths[q],
		       vox_origin, vox_size, v_coords[0], h_pixels[p+q],
		       cphi, sphi, nx, ny, nz, "pq");
	  }
	  if (p+cnt < nh_pixels) {
	    test_voxel(i, j, 0, a, 0, p+cnt, 0.0, 0.0,
		       vox_origin, vox_size, v_coords[0], h_pixels[p+cnt],
		       cphi, sphi, nx, ny, nz, "p+");
	  }
#endif // TESTBP
	}
      }
    }
  }
}
