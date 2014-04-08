
// Ugly copy of projection algorithm, could simplify to 2D
// Todo - how to avoid a copy?
#include <map>
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#ifdef MATLAB_MEX_FILE
#  include "mex_types.hpp"
#else
#  include "base_types.hpp"
#endif // mex
#include "instruments.hpp"

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

#define PRECISION real(0.00000001) /* for calculating rays intersecting voxels*/

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

void CCPi::parallel_beam::map_2Dprojection(const real start[], const real end[],
					   const real b_x, const real b_y,
					   const real b_z, const real d_x,
					   const real d_y, const real d_z,
					   const int im_size_x,
					   const int im_size_y,
					   const int im_size_z,
					   const sl_int z_offset,
					   projection_map &map)
{

  int N_x, N_y, N_z, N_p/*, im_size_x, im_size_y, im_size_z*/;
  real /*b_x, b_y, b_z, d_x, d_y, d_z,*/ d_conv;
    real p1_x, p1_y, p1_z, p2_x, p2_y, p2_z;

    int x_defined, y_defined, z_defined;
    sl_int i=0,j=0,k=0;

    real alpha_x_min, alpha_y_min, alpha_z_min, alpha_x_max, alpha_y_max,
	alpha_z_max, alpha_min, alpha_max, alpha_x, alpha_y, alpha_z, alpha_c;
    real alpha_x_u = 0.0, alpha_y_u = 0.0, alpha_z_u = 0.0;
    real l_ij;
    int i_min, j_min, k_min, i_max, j_max, k_max, n_count, i_u, j_u, k_u;
    sl_int i_step, j_step, k_step;

    sl_int ray_index;

    p1_x = start[0];
    p1_y = start[1];
    p1_z = start[2];
    p2_x = end[0];
    p2_y = end[1];
    p2_z = end[2];

    //im_size_x = options->im_size_x;
    ///im_size_y = options->im_size_y;
    //im_size_z = options->im_size_z;

	N_x=im_size_x+1;
	N_y=im_size_y+1;
	N_z=im_size_z+1;

	/* d: voxel size */
	//d_x = options->d_x;
	//d_y = options->d_y;
	//d_z = options->d_z;

	/* b: grid offset from origin */
	//b_x = options->b_x;
	//b_y = options->b_y;
	//b_z = options->b_z;



    /* use total lengh=alpha_max-alpha_min instead, to get everage, not sum. */
    /* moving back to original d_conv*/
    d_conv=sqrt( (p1_x-p2_x)*(p1_x-p2_x) + (p1_y-p2_y)*(p1_y-p2_y) + (p1_z-p2_z)*(p1_z-p2_z));


    x_defined =  !(equal_to_precision(p1_x,p2_x,PRECISION));
    y_defined =  !(equal_to_precision(p1_y,p2_y,PRECISION));
    z_defined =  !(equal_to_precision(p1_z,p2_z,PRECISION));

    if( !x_defined && !y_defined && !z_defined)
	return;

    if (x_defined) {
	alpha_x_min=std::min(alpha_fn(0, p1_x, p2_x, b_x, d_x), alpha_fn(N_x-1, p1_x, p2_x, b_x, d_x));
	alpha_x_max=std::max(alpha_fn(0, p1_x, p2_x, b_x, d_x), alpha_fn(N_x-1, p1_x, p2_x, b_x, d_x));
    }
    else {
	alpha_x_min=-2;
	alpha_x_max=2;
	i=(int) floor_j( phix(real(0.0), p1_x, p2_x, b_x, d_x));
	if ( i < 0 || i >= im_size_x)
	    return;
	alpha_x=2;
	i_min = 1;
	i_max = 0;
    }

    if(y_defined) {
	alpha_y_min=std::min(alpha_fn(0, p1_y, p2_y, b_y, d_y), alpha_fn(N_y-1, p1_y, p2_y, b_y, d_y));
	alpha_y_max=std::max(alpha_fn(0, p1_y, p2_y, b_y, d_y), alpha_fn(N_y-1, p1_y, p2_y, b_y, d_y));
    }
    else {
	alpha_y_min=-2;
	alpha_y_max=2;
	j=(int) floor_j( phix(real(0.0), p1_y, p2_y, b_y, d_y));
	if ( j < 0 || j >= im_size_y)
	    return;
	alpha_y=2;
	j_min = 1;
	j_max = 0;
    }


    if(z_defined) {
	alpha_z_min=std::min(alpha_fn(0, p1_z, p2_z, b_z, d_z), alpha_fn(N_z-1, p1_z, p2_z, b_z, d_z));
	alpha_z_max=std::max(alpha_fn(0, p1_z, p2_z, b_z, d_z), alpha_fn(N_z-1, p1_z, p2_z, b_z, d_z));
    }
    else {
	alpha_z_min=-2;
	alpha_z_max=2;
	k=(int) floor_j( phix(real(0.0), p1_z, p2_z, b_z, d_z));
	if ( k < 0 || k >= im_size_z)
	    return;
	alpha_z=2;
	k_min = 1;
	k_max = 0;
    }

    alpha_min=std::max(real(0.0),
		       max3_dbl(alpha_x_min, alpha_y_min, alpha_z_min));
    alpha_max=std::min(real(1.0),
		       min3_dbl(alpha_x_max, alpha_y_max, alpha_z_max));

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



	if (z_defined && p1_z < p2_z) {
	    if (equal_to_precision(alpha_min,alpha_z_min,PRECISION)==1)
		k_min=1;
	    else
		k_min = (int) ceil_j(phix(alpha_min, p1_z, p2_z, b_z, d_z));


	    if (equal_to_precision(alpha_max, alpha_z_max,PRECISION)==1)
		k_max = N_z - 1;
	    else
		k_max = (int) floor_j( phix(alpha_max, p1_z, p2_z, b_z, d_z));

	    alpha_z=alpha_fn(k_min, p1_z, p2_z, b_z, d_z);
	}

	else if (z_defined) {

	    if (equal_to_precision(alpha_min,alpha_z_min,PRECISION)==1)
		k_max=N_z-2;
	    else
		k_max = (int) floor_j(phix(alpha_min, p1_z, p2_z, b_z, d_z));


	    if (equal_to_precision(alpha_max, alpha_z_max, PRECISION)==1)
		k_min = 0;
	    else
		k_min = (int) ceil_j( phix(alpha_max, p1_z, p2_z, b_z, d_z));

	    alpha_z=alpha_fn(k_max, p1_z, p2_z, b_z, d_z);
	}


	N_p=(i_max - i_min +1) + (j_max - j_min + 1) + (k_max - k_min + 1);

	if (x_defined) {
	    i=(int) floor_j( phix( (min3_dbl(alpha_x, alpha_y, alpha_z) + alpha_min)/2, p1_x, p2_x, b_x, d_x) );
	alpha_x_u = d_x/std::abs(p2_x-p1_x);
	}

	if (y_defined) {
	    j=(int) floor_j( phix( (min3_dbl(alpha_x, alpha_y, alpha_z) + alpha_min)/2, p1_y, p2_y, b_y, d_y) );
	alpha_y_u = d_y/std::abs(p2_y-p1_y);
	}
	if (z_defined) {
	    k=(int) floor_j( phix( (min3_dbl(alpha_x, alpha_y, alpha_z) + alpha_min)/2, p1_z, p2_z, b_z, d_z) );
	alpha_z_u = d_z/std::abs(p2_z-p1_z);
	}

	if (p1_x < p2_x)
	    i_u=1;
	else
	    i_u=-1;

	if (p1_y < p2_y)
	    j_u=1;
	else
	    j_u=-1;

	if (p1_z < p2_z)
	    k_u=1;
	else
	    k_u=-1;


	alpha_c=alpha_min;
	ray_index = (k+z_offset)*im_size_y*im_size_x + j*im_size_x + i;
	i_step = i_u;
	j_step = j_u * im_size_x;
	k_step = k_u * im_size_y * im_size_x;
	map_index index;
	real data = 0.0;

	for (n_count=1; n_count<N_p+1;n_count++) {


	    /* x smallest*/
	    if (x_defined && alpha_x <= alpha_y && alpha_x <= alpha_z) {
		/* ray intersects pixel(i,j) with length l_ij */

	      data = (alpha_x-alpha_c)*d_conv;
	      index.x = j;
	      index.y = i;
	      map.insert(std::make_pair(index, data));

		if( y_defined && alpha_x == alpha_y) {
		    j += j_u;
		    ray_index += j_step;
		    n_count++;
		    alpha_y += alpha_y_u;
		}

		if( z_defined && alpha_x == alpha_z) {
		    k += k_u;
		    ray_index += k_step;
		    n_count++;
		    alpha_z += alpha_z_u;
		}

		i += i_u;
		ray_index += i_step;
		alpha_c=alpha_x;
		alpha_x += alpha_x_u;
	    }

	    /* y smallest*/
	    else if (y_defined && alpha_y <= alpha_z) {
		/* ray intersects pixel(i,j) with length l_ij */

	      data = (alpha_y-alpha_c)*d_conv;
	      index.x = j;
	      index.y = i;
	      map.insert(std::make_pair(index, data));

		if( z_defined && alpha_y == alpha_z) {
		    k += k_u;
		    ray_index += k_step;
		    n_count++;
		    alpha_z += alpha_z_u;
		}

		j=j+j_u;
		ray_index += j_step;
		alpha_c=alpha_y;
		alpha_y += alpha_y_u;
	    }

	    /* z smallest*/
	    else if (z_defined) {
		/* ray intersects pixel(i,j) with length l_ij */

	      data = (alpha_z-alpha_c)*d_conv;
	      index.x = j;
	      index.y = i;
	      map.insert(std::make_pair(index, data));

		k += k_u;
		ray_index += k_step;
		alpha_c=alpha_z;
		alpha_z += alpha_z_u;
	    }


	    /* did we loop too far? */
	    if( i < 0 || j < 0 || k < 0 || i >= im_size_x || j >= im_size_y || k >= im_size_z)
		/* artificially end loop  */
		N_p = n_count - 1;



	} /* end of for loop though N_p */

	/* in case we're ending inside grid, finish off last voxel */
	if( (alpha_max - alpha_c) > PRECISION) {
	    /* this is the last step so don't need to worry about incrementing i or j*/
	    l_ij=(alpha_max-alpha_c)*d_conv;

	    data = l_ij;
	    index.x = j;
	    index.y = i;
	    map.insert(std::make_pair(index, data));
	}

    } /* of alpha_min < alpha_max */

    return;
}

void CCPi::parallel_beam::setup_2D_matrix(const real det_y[], const real phi[],
					  const int n_angles,
					  const int n_rays_z,
					  const int n_rays_y,
					  const real grid_offset[3],
					  const real voxel_size[3],
					  const int nx_voxels,
					  const int ny_voxels,
					  const int nz_voxels)
{
  sl_int curr_angle, curr_ray_y;
  real cos_curr_angle, sin_curr_angle;
  real start[3], end[3];

  // set detector z to 2* the yz limits of the voxels, so it misses
  // longest voxel dim should be sqrt(3), so 2 should be safe
  real det_x = real(2.0) * std::max(std::abs(grid_offset[0]),
				    std::max(std::abs(grid_offset[1]),
					     std::abs(grid_offset[2])));

  end[2] = 0.000001;
  start[2] = end[2];

  // generate forward projections.
  real **forward_data = new real *[n_rays_y * n_angles];
  int **forward_x = new int *[n_rays_y * n_angles];
  int **forward_y = new int *[n_rays_y * n_angles];
  int *forward_sizes = new int[n_rays_y * n_angles];

  sl_int total = 0;
  // Todo parallelize?
  int offset = 0;
  for (curr_angle = 0; curr_angle < n_angles; curr_angle++) {
    /* rotate source and detector positions by current angle */
    cos_curr_angle = std::cos(phi[curr_angle]);
    sin_curr_angle = std::sin(phi[curr_angle]);

    /* loop over y values on detector */
    for (curr_ray_y = 0; curr_ray_y < n_rays_y; curr_ray_y++) {
      end[0] = cos_curr_angle * det_x - sin_curr_angle * det_y[curr_ray_y];
      end[1] = sin_curr_angle * det_x + cos_curr_angle * det_y[curr_ray_y];
      start[0] = end[0] - real(3.0) * cos_curr_angle * det_x;
      start[1] = end[1] - real(3.0) * sin_curr_angle * det_x;

      /* loop over z values on detector */

      projection_map map;
      map_2Dprojection(start, end, grid_offset[0], grid_offset[1],
		       grid_offset[2], voxel_size[0], voxel_size[1],
		       voxel_size[2], nx_voxels, ny_voxels, nz_voxels, 0,
		       map);
      // Do stuff with map.
      int sz = (int)map.size();
      forward_sizes[offset] = sz;
      if (sz > 0) {
	total += sz;
	forward_data[offset] = new real[sz];
	forward_x[offset] = new int[sz];
	forward_y[offset] = new int[sz];
	int i = 0;
	for (projection_map::const_iterator ptr = map.begin();
	     ptr != map.end(); ++ptr) {
	  forward_data[offset][i] = ptr->second;
	  forward_x[offset][i] = ptr->first.y;
	  forward_y[offset][i] = ptr->first.x;
	  i++;
	}
      } else {
	forward_data[offset] = 0;
	forward_x[offset] = 0;
	forward_y[offset] = 0;
      }
      offset++;
    }
  }

  // forward matrix - the pixels produced have to be offset by z
  // MKL 0 indexed coordinate format, CSR is a bit tricky as needs all rows
  // even though we only want a 2D subset (its because of z in angles,z,y)
  matrix_size = total;
  forward_matrix = new real[total];
  forward_cols = new sl_int[total];
  forward_rows = new sl_int[total];
  sl_int column = 0;
  sl_int row = 0;
  sl_int index = 0;
  offset = 0;
  for (curr_angle = 0; curr_angle < n_angles; curr_angle++) {
    for (curr_ray_y = 0; curr_ray_y < n_rays_y; curr_ray_y++) {
      row = curr_angle * n_rays_y * n_rays_z + curr_ray_y;
      if (forward_sizes[offset] > 0) {
	for (int i = 0; i < forward_sizes[offset]; i++) {
	  column = forward_y[offset][i] * nx_voxels + forward_x[offset][i];
	  forward_matrix[index] = forward_data[offset][i];
	  forward_cols[index] = column;
	  forward_rows[index] = row;
	  index++;
	}
	
      }
      offset++;
    }
  }

  // produce backward projections from forward data - since we have packed
  // y,x memory order we can use CSR format
  backward_matrix = new real[total];
  backward_cols = new sl_int[total];
  backward_rowb = new sl_int[nx_voxels * ny_voxels];
  backward_rowe = new sl_int[nx_voxels * ny_voxels];

  sl_int count = 0;
  column = 0;
  row = 0;
  index = 0;
  offset = 0;
  // store last position to reduce search loop
  std::vector<int> position(n_angles * n_rays_y);
  for (int i = 0; i < n_angles * n_rays_y; i++)
    position[i] = 0;
  for (sl_int y = 0; y < ny_voxels; y++) {
    std::vector<projection_map> mapping(nx_voxels);
    int p = 0;
    map_index m;
    for (curr_angle = 0; curr_angle < n_angles; curr_angle++) {
      for (curr_ray_y = 0; curr_ray_y < n_rays_y; curr_ray_y++) {
	for (int i = position[p]; i < forward_sizes[p]; i++) {
	  if (forward_y[p][i] == y) {
	    m.x = curr_angle;
	    m.y = curr_ray_y;
	    mapping[forward_x[p][i]].insert(std::make_pair(m,
							   forward_data[p][i]));
	  } else if (forward_y[p][i] > y) {
	    position[p] = i;
	    break;
	  }
	}
	p++;
      }
    }

    for (sl_int x = 0; x < nx_voxels; x++) {
      int sz = (int)mapping[x].size();
      backward_rowb[offset] = index;
      if (sz > 0) {
	count += sz;
	for (projection_map::const_iterator ptr = mapping[x].begin();
	     ptr != mapping[x].end(); ++ptr) {
	  column = sl_int(ptr->first.x) * sl_int(n_rays_y) * sl_int(n_rays_z)
	    + sl_int(ptr->first.y);
	  backward_matrix[index] = ptr->second;
	  backward_cols[index] = column;
	  index++;
	}
      }
      backward_rowe[offset] = index;
      offset++;
    }
  }
  for (int i = 0; i < n_angles * n_rays_y; i++) {
    if (forward_sizes[i] > 0) {
      delete [] forward_y[i];
      delete [] forward_x[i];
      delete [] forward_data[i];
    }
  }
  delete [] forward_sizes;
  delete [] forward_y;
  delete [] forward_x;
  delete [] forward_data;
#ifdef USE_TIMER
  std::cout << "Total matrix data is " << total << ' ' << count << '\n';
#endif // USE_TIMER
}
