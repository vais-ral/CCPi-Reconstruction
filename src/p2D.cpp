
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

extern void test_voxel(const int i, const int j, const int k,
		       const int a, const int v, const int h, const double d,
		       const double l, const double voxel_origin[3],
		       const double voxel_size[3], const double v_pos,
		       const double h_pos, const double cphi,
		       const double sphi, const int nx, const int ny,
		       const int nz, const char message[]);

/* jacobs_ray_3d
   % void jacobs_ray_3d(int im_size, real *start, real *end, int *ray_index, real *ray_data, int *n_entries)
   %
   % implementation of algorithm by Jacobs et.al (modification of Siddon's)
   %
   % int im_size: size of grid
   % real start[3]: x,y,z coordinates for starting point of ray
   % real end[3]
   % int ray_index[n_entries]: stores voxel numbers hit by ray
   % real ray_data[n_entries]: stores weights corresponding to relative length of ray in a given voxel
   % int n_entries: counts number of voxels hit by ray
   %
   %
   % takes coordinates of points relative origin in centre of grid, and
   % translates to correspond to algorithm
   %
   %
   % david szotten
   % august 2005

$Revision: 1.1.1.11 $
$Date: 2008/09/08 13:20:38 $
*/

/* 18/08/2011 WT - adapting this code to perform back projection (ie At*b) in single step */

/* 09/09/2011 this version operates on single precision ray and volume data, but performs internal calculations in double precision */

#define PRECISION 0.00000001 /* for calculating rays intersecting voxels*/

static inline real alpha_fn(const int n, const real p1, const real p2,
			    const real b, const real d)
{
    return ( (b+n*d) - p1)/(p2-p1);
}

static inline real p(const real alpha, const real p1, const real p2)
{
    return p1+alpha*(p2-p1);
}

static inline real phix(const real alpha, const real p1, const real p2,
			const real b, const real d)
{
    return ( p(alpha, p1, p2)-b)/d;
}

static inline int equal_to_precision(const real x, const real y,
				     const real prec)
{
    return std::abs(x-y) < prec;
}

static inline real min3_dbl(const real a, const real b, const real c)
{
  return a < b ? std::min(a,c) : std::min(b,c);
}

static inline real max3_dbl(const real a, const real b, const real c)
{
  return a > b ? std::max(a,c) : std::max(b,c);
}

static inline real ceil_j(const real arg)
{
  return arg == (int)arg ? arg+1 : std::ceil( arg );
}

static inline real floor_j(const real arg)
{
  return std::floor( arg );
}

void CCPi::parallel_beam::f2D_parallel(const real start[], const real end[],
				       const real b_x, const real b_y,
				       const real d_x, const real d_y,
				       const int im_size_x,
				       const int im_size_y,
				       const int im_size_z,
				       voxel_type *const voxels,
				       pixel_type *pixels,
				       const int pixels_per_voxel)
{

  int N_x, N_y, N_p;
  real d_conv;
  real p1_x, p1_y, p2_x, p2_y;

  int x_defined, y_defined;
  long i=0,j=0;

  real alpha_x_min, alpha_y_min, alpha_x_max, alpha_y_max,
    alpha_min, alpha_max, alpha_x, alpha_y, alpha_c;
  real alpha_x_u = 0.0, alpha_y_u = 0.0;
  real l_ij;
  int i_min, j_min, i_max, j_max, n_count, i_u, j_u;
  long i_step, j_step;

  long ray_index;

  p1_x = start[0];
  p1_y = start[1];
  p2_x = end[0];
  p2_y = end[1];

  N_x=im_size_x+1;
  N_y=im_size_y+1;

  /* use total lengh=alpha_max-alpha_min instead, to get everage, not sum. */
  /* moving back to original d_conv*/
  d_conv=sqrt( (p1_x-p2_x)*(p1_x-p2_x) + (p1_y-p2_y)*(p1_y-p2_y));

  x_defined =  !(equal_to_precision(p1_x,p2_x,PRECISION));
  y_defined =  !(equal_to_precision(p1_y,p2_y,PRECISION));

  if( !x_defined && !y_defined)
    return;

  if (x_defined) {
    alpha_x_min=std::min(alpha_fn(0, p1_x, p2_x, b_x, d_x), alpha_fn(N_x-1, p1_x, p2_x, b_x, d_x));
    alpha_x_max=std::max(alpha_fn(0, p1_x, p2_x, b_x, d_x), alpha_fn(N_x-1, p1_x, p2_x, b_x, d_x));
  } else {
    alpha_x_min=-2;
    alpha_x_max=2;
    i=(int) floor_j( phix(0.0, p1_x, p2_x, b_x, d_x));
    if ( i < 0 || i >= im_size_x)
      return;
    alpha_x=2;
    i_min = 1;
    i_max = 0;
  }

  if(y_defined) {
    alpha_y_min=std::min(alpha_fn(0, p1_y, p2_y, b_y, d_y), alpha_fn(N_y-1, p1_y, p2_y, b_y, d_y));
    alpha_y_max=std::max(alpha_fn(0, p1_y, p2_y, b_y, d_y), alpha_fn(N_y-1, p1_y, p2_y, b_y, d_y));
  } else {
    alpha_y_min=-2;
    alpha_y_max=2;
    j=(int) floor_j( phix(0.0, p1_y, p2_y, b_y, d_y));
    if ( j < 0 || j >= im_size_y)
      return;
    alpha_y=2;
    j_min = 1;
    j_max = 0;
  }

  alpha_min=std::max(0.0, std::max(alpha_x_min, alpha_y_min));
  alpha_max=std::min(1.0, std::min(alpha_x_max, alpha_y_max));

  /* if ray intersects voxel grid */
  if (alpha_min < alpha_max) {

    if (x_defined && p1_x < p2_x) {
      if (equal_to_precision(alpha_min,alpha_x_min,PRECISION)==1)
	i_min=1;
      else
	i_min = (int) ceil_j(phix(alpha_min, p1_x, p2_x, b_x, d_x));

      if (equal_to_precision(alpha_max,alpha_x_max,PRECISION)==1)
	i_max = N_x - 1;
      else
	i_max = (int) floor_j( phix(alpha_max, p1_x, p2_x, b_x, d_x));
      
      alpha_x=alpha_fn(i_min, p1_x, p2_x, b_x, d_x);
    }

    else if (x_defined) {
      if (equal_to_precision(alpha_min,alpha_x_min,PRECISION)==1)
	i_max=N_x-2;
      else
	i_max = (int) floor_j(phix(alpha_min, p1_x, p2_x, b_x, d_x));

      if (equal_to_precision(alpha_max,alpha_x_max,PRECISION)==1)
	i_min = 0;
      else
	i_min = (int) ceil_j( phix(alpha_max, p1_x, p2_x, b_x, d_x));
      
      alpha_x=alpha_fn(i_max, p1_x, p2_x, b_x, d_x);
    }

    if (y_defined && p1_y < p2_y) {
      if (equal_to_precision(alpha_min,alpha_y_min,PRECISION)==1)
	j_min=1;
      else
	j_min = (int) ceil_j(phix(alpha_min, p1_y, p2_y, b_y, d_y));

      if (equal_to_precision(alpha_max, alpha_y_max,PRECISION)==1)
	j_max = N_y - 1;
      else
	j_max = (int) floor_j( phix(alpha_max, p1_y, p2_y, b_y, d_y));

      alpha_y=alpha_fn(j_min, p1_y, p2_y, b_y, d_y);
    }

    else if (y_defined) {

      if (equal_to_precision(alpha_min,alpha_y_min,PRECISION)==1)
	j_max=N_y-2;
      else
	j_max = (int) floor_j(phix(alpha_min, p1_y, p2_y, b_y, d_y));
      
      if (equal_to_precision(alpha_max, alpha_y_max, PRECISION)==1)
	j_min = 0;
      else
	j_min = (int) ceil_j( phix(alpha_max, p1_y, p2_y, b_y, d_y));
      
      alpha_y=alpha_fn(j_max, p1_y, p2_y, b_y, d_y);
    }

    N_p=(i_max - i_min +1) + (j_max - j_min + 1);

    if (x_defined) {
      i=(int) floor_j( phix( (std::min(alpha_x, alpha_y) + alpha_min)/2, p1_x, p2_x, b_x, d_x) );
      alpha_x_u = d_x/std::abs(p2_x-p1_x);
    }

    if (y_defined) {
      j=(int) floor_j( phix( (std::min(alpha_x, alpha_y) + alpha_min)/2, p1_y, p2_y, b_y, d_y) );
      alpha_y_u = d_y/std::abs(p2_y-p1_y);
    }

    if (p1_x < p2_x)
      i_u=1;
    else
      i_u=-1;

    if (p1_y < p2_y)
      j_u=1;
    else
      j_u=-1;

    alpha_c=alpha_min;
    ray_index = i*im_size_y*im_size_z + j*im_size_z;
    i_step = i_u * im_size_y * im_size_z;
    j_step = j_u * im_size_z;
    real data = 0.0;

    for (n_count=1; n_count<N_p+1;n_count++) {

      /* x smallest*/
      if (x_defined && alpha_x <= alpha_y) {
	/* ray intersects pixel(i,j) with length l_ij */

	data = (alpha_x-alpha_c)*d_conv;
	/*
	// Some of this can also be turned into a vector op with proper align
	switch (pixels_per_voxel) {
	case 1:
	  for (int k = 0; k < im_size_z; k++) {
	    pixels[k] += data * voxels[ray_index + k];
	  }
	  break;
	case 2:
	  int k_range = 0;
	  for (int k = 0; k < im_size_z; k++) {
	    pixels[k_range + 0] += data * voxels[ray_index + k];
	    pixels[k_range + 1] += data * voxels[ray_index + k];
	    k_range += 2;
	  }
	  break;
	case 4:
	  int k_range = 0;
	  for (int k = 0; k < im_size_z; k++) {
	    for (int kp = 0; kp < 4; kp++)
	      pixels[k_range + kp] += data * voxels[ray_index + k];
	    k_range += 4;
	  }
	  break;
	default:
	  int k_range = 0;
	  for (int k = 0; k < im_size_z; k++) {
	    for (int kp = 0; kp < pixels_per_voxel; kp++)
	      pixels[k_range + kp] += data * voxels[ray_index + k];
	    k_range += pixels_per_voxel;
	  }
	  break;
	}
	*/

	int k_range = 0;
	for (int k = 0; k < im_size_z; k++) {
	  for (int kp = 0; kp < pixels_per_voxel; kp++)
	    pixels[k_range + kp] += data * voxels[ray_index + k];
	  k_range += pixels_per_voxel;
	}

	if( y_defined && alpha_x == alpha_y) {
	  j += j_u;
	  ray_index += j_step;
	  n_count++;
	  alpha_y += alpha_y_u;
	}
	i += i_u;
	ray_index += i_step;
	alpha_c=alpha_x;
	alpha_x += alpha_x_u;
      }
      
      /* y smallest*/
      else if (y_defined) {
	/* ray intersects pixel(i,j) with length l_ij */

	data = (alpha_y-alpha_c)*d_conv;
	int k_range = 0;
	for (int k = 0; k < im_size_z; k++) {
	  for (int kp = 0; kp < pixels_per_voxel; kp++)
	    pixels[k_range + kp] += data * voxels[ray_index + k];
	  k_range += pixels_per_voxel;
	}
	j=j+j_u;
	ray_index += j_step;
	alpha_c=alpha_y;
	alpha_y += alpha_y_u;
      }
      
      /* did we loop too far? */
      if( i < 0 || j < 0 || i >= im_size_x || j >= im_size_y)
	/* artificially end loop  */
	N_p = n_count - 1;
    } /* end of for loop though N_p */
    
    /* in case we're ending inside grid, finish off last voxel */
    if( (alpha_max - alpha_c) > PRECISION) {
      /* this is the last step so don't need to worry about incrementing i or j*/
      l_ij=(alpha_max-alpha_c)*d_conv;
      
      data = l_ij;
      int k_range = 0;
      for (int k = 0; k < im_size_z; k++) {
	for (int kp = 0; kp < pixels_per_voxel; kp++)
	  pixels[k_range + kp] += data * voxels[ray_index + k];
	k_range += pixels_per_voxel;
      }
    }
    
  } /* of alpha_min < alpha_max */
  
  return;
}

void CCPi::parallel_beam::f2D(const real det_y[], const real phi[],
			      const int n_angles,
			      const int n_rays_y,
			      const int n_rays_z,
			      const real grid_offset[3],
			      const real voxel_size[3],
			      const int nx_voxels,
			      const int ny_voxels,
			      const int nz_voxels,
			      pixel_type ray_data[],
			      voxel_type *const vol_data) const
{
  long curr_angle, curr_ray_y;
  real cos_curr_angle, sin_curr_angle;
  real start[3], end[3];

  // set detector z to 2* the yz limits of the voxels, so it misses
  // longest voxel dim should be sqrt(3), so 2 should be safe
  real det_x = 2.0 * std::max(std::abs(grid_offset[0]),
			      std::max(std::abs(grid_offset[1]),
				       std::abs(grid_offset[2])));

  const int pixels_per_voxel = n_rays_z / nz_voxels;
#pragma omp parallel for shared(det_y, ray_data, phi) private(curr_angle, curr_ray_y, start, end, cos_curr_angle, sin_curr_angle), firstprivate(det_x, pixels_per_voxel) schedule(dynamic)
  for (curr_angle = 0; curr_angle < n_angles; curr_angle++) {
    //end[2] = 0.000001;
    //start[2] = end[2];
    /* rotate source and detector positions by current angle */
    cos_curr_angle = std::cos(phi[curr_angle]);
    sin_curr_angle = std::sin(phi[curr_angle]);

    long a_offset = curr_angle * n_rays_y * n_rays_z;
    /* loop over y values on detector */
    for (curr_ray_y = 0; curr_ray_y < n_rays_y; curr_ray_y++) {
      end[0] = cos_curr_angle * det_x - sin_curr_angle * det_y[curr_ray_y];
      end[1] = sin_curr_angle * det_x + cos_curr_angle * det_y[curr_ray_y];
      start[0] = end[0] - 3.0 * cos_curr_angle * det_x;
      start[1] = end[1] - 3.0 * sin_curr_angle * det_x;

      long y_offset = a_offset + curr_ray_y * n_rays_z;

      f2D_parallel(start, end, grid_offset[0], grid_offset[1],
		   voxel_size[0], voxel_size[1], nx_voxels, ny_voxels,
		   nz_voxels, vol_data, &ray_data[y_offset], pixels_per_voxel);
    }
  }
}

// Todo - reorder i,j; think about scaling to 1.0 ideas, best loop order...
// Loop blocking to improve cache reuse?
void CCPi::parallel_beam::my_back_project(const real h_pixels[],
					  const real v_pixels[],
					  const real angles[],
					  pixel_type pixels[],
					  voxel_type *const voxels,
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
  long nyz = long(ny) * long(nz);
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
  // Todo - pre-rotated array of voxels? need for each angle, x, y though
  // x' = x cos() + y sin()
  // y' = -x sin() + y cos() - not used in current form
  std::vector<real> vox_lengths(n_angles);
  for (int i = 0; i < n_angles; i++) {
    if (std::abs(cangles[i]) > std::abs(sangles[i]))
      vox_lengths[i] = vox_size[0] / std::abs(cangles[i]);
    else
      vox_lengths[i] = vox_size[0] / std::abs(sangles[i]);
  }
  // Todo can precalc 2D i0/inh1, yc0, xh_coords here
  // its probably p that is the tricky one though
  // put end value to stop while loop without over-run
  std::vector<real> v_coords(nv_pixels + 1);
  for (int i = 0; i < nv_pixels; i++)
    v_coords[i] = v_pixels[i];
  v_coords[nv_pixels] = z_coords[nz] + 1.0;
  std::vector<real> xx_coords(nx + 1);
  std::vector<real> yy_coords(ny + 1);
  if (v_pixels[0] < z_coords[0] or v_pixels[pixels_per_voxel] < z_coords[1])
    std::cerr << "Oops - bad back projection\n";
  // Todo - what is best order of a/j loops and which to parallelise?
  for (int a = 0; a < n_angles; a++) {
    long pix_av = a * long(nv_pixels) * long(nh_pixels);
    real cphi = cangles[a];
    real sphi = sangles[a];
    real L = vox_lengths[a];
    //for (int i = 0; i < nx; i++)
    //xh_coords[i] = - x_coords[i] * sphi * inv_pixel_step;
    for (int i = 0; i <= nx; i++)
      xx_coords[i] = - x0_coords[i] * sphi;
    for (int i = 0; i <= ny; i++)
      yy_coords[i] = y0_coords[i] * cphi;
    real y_00 = xx_coords[0] + yy_coords[0];
    real y_01 = xx_coords[0] + yy_coords[1];
    real y_10 = xx_coords[1] + yy_coords[0];
    real y_11 = xx_coords[1] + yy_coords[1];
    real y_top = 0.0;
    real y_l1;
    real y_l2;
    real y_bot = 0.0;
    if (sphi >= 0.0) {
      // <= 180, scanning from 0 to nx decreases p
      // orientation of voxels is the same for all so can get rotated shape
      // from any single one - use 0.0
      int jm;
      if (cphi >= 0.0) { // 0-90
	jm = 1;
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
	jm = 0;
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
      //const int npixels = int(y_bot * inv_pixel_step) + 2;
      // Its always scaled by L
      real inv_y_l1 = L / y_l1;
      //#pragma omp parallel for shared(xx_coords, yy_coords, h_pixels, pixels) firstprivate(inv_y_l1, pixels_per_voxel, nx, ny, nz, nh_pixels, inv_pixel_step, pix_av, y_l1, y_l2, y_top, y_bot, nyz, cphi, sphi, L) schedule(dynamic)
      for (int j = 0; j < ny; j++) {
#ifdef TESTBP
	std::vector<real> lengths(nh_pixels);
#endif // TESTBP
	real yj = yy_coords[j + jm];
	long vox_zy = j * long(nz);
	int p = int((xx_coords[0] + yj - h_pixels[0]) * inv_pixel_step);
	if (p >= nh_pixels)
	  p = nh_pixels - 1;
	else if (p < 0)
	  continue;
	for (int i = 0; i < nx; i++) {
	  // Todo ? ytop decreases by y_l1 for each i step
	  long vox_xy = vox_zy + i * nyz;
	  y_top = xx_coords[i] + yj;
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
	      int k_range = pix_av + q * long(nv_pixels);
	      for (int k = 0; k < nz; k++) {
		for (int v = 0; v < pixels_per_voxel; v++)
		  voxels[vox_xy + k] += pixels[k_range + v] * len;
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
	}
      }
    } else {
      // 0-nx increases p
      int jm;
      if (cphi < 0.0) { // 180-270
	jm = 1;
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
	jm = 0;
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
      //#pragma omp parallel for shared(xx_coords, yy_coords, h_pixels, pixels) firstprivate(inv_y_l1, pixels_per_voxel, nx, ny, nz, nh_pixels, inv_pixel_step, pix_av, y_l1, y_l2, y_top, y_bot, nyz, cphi, sphi, L) schedule(dynamic)
      for (int j = 0; j < ny; j++) {
#ifdef TESTBP
	std::vector<real> lengths(nh_pixels);
#endif // TESTBP
	real yj = yy_coords[j + jm];
	long vox_zy = j * long(nz);
	int p = int((xx_coords[0] + yj - h_pixels[0]) * inv_pixel_step);
	if (p < 0)
	  p = 0;
	else if (p >= nh_pixels)
	  continue;
	for (int i = 0; i < nx; i++) {
	  long vox_xy = vox_zy + i * nyz;
	  // Todo ? ybot increases by y_l1 for each i step
	  y_bot = xx_coords[i] + yj;
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
	      int k_range = pix_av + q * long(nv_pixels);
	      for (int k = 0; k < nz; k++) {
		for (int v = 0; v < pixels_per_voxel; v++)
		  voxels[vox_xy + k] += pixels[k_range + v] * len;
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
