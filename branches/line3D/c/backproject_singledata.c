/* jacobs_ray_3d
   % void jacobs_ray_3d(int im_size, double *start, double *end, int *ray_index, double *ray_data, int *n_entries)
   %
   % implementation of algorithm by Jacobs et.al (modification of Siddon's)
   %
   % int im_size: size of grid
   % double start[3]: x,y,z coordinates for starting point of ray
   % double end[3]
   % int ray_index[n_entries]: stores voxel numbers hit by ray
   % double ray_data[n_entries]: stores weights corresponding to relative length of ray in a given voxel
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

#include "jacobs_rays.h"
#include "omp.h"

static inline double alpha_fn(const int n, const double p1, const double p2,
			      const double b, const double d)
{
    return ( (b+n*d) - p1)/(p2-p1);
}

static inline double p(const double alpha, const double p1, const double p2)
{
    return p1+alpha*(p2-p1);
}

static inline double phi(const double alpha, const double p1, const double p2,
			 const double b, const double d)
{
    return ( p(alpha, p1, p2)-b)/d;
}

static inline int equal_to_precision(const double x, const double y,
				     const double prec)
{
    return fabs(x-y) < prec;
}

static inline double min_dbl(const double a, const double b)
{
  return a < b ? a : b;
}

static inline double max_dbl(const double a, const double b)
{
  return a > b ? a : b;
}

static inline double min3_dbl(const double a, const double b, const double c)
{
  return a < b ? min_dbl(a,c) : min_dbl(b,c);
}

static inline double max3_dbl(const double a, const double b, const double c)
{
  return a > b ? max_dbl(a,c) : max_dbl(b,c);
}

static inline double ceil_j(const double arg)
{
  return arg == (int)arg ? arg+1 : ceil( arg );
}

static inline double floor_j(const double arg)
{
  return floor( arg );
}

void backproject_singledata(const double start[], const double end[],
			    const double ray_data, REAL vol_data[],
			    const struct jacobs_options *options,
			    const long z_offset)
{
    
    int N_x, N_y, N_z, N_p, im_size_x, im_size_y, im_size_z;
    double b_x, b_y, b_z, d_x, d_y, d_z, d_conv;
    double p1_x, p1_y, p1_z, p2_x, p2_y, p2_z;
    
    int x_defined, y_defined, z_defined;
    long i,j,k;
    
    double alpha_x_min, alpha_y_min, alpha_z_min, alpha_x_max, alpha_y_max, 
	alpha_z_max, alpha_min, alpha_max, alpha_x, alpha_y, alpha_z, alpha_c;
    double alpha_x_u, alpha_y_u, alpha_z_u;
    double l_ij;
    int i_min, j_min, k_min, i_max, j_max, k_max, n_count, i_u, j_u, k_u;
    long i_step, j_step, k_step;
    
    long ray_index;
	
    p1_x = start[0];
    p1_y = start[1];
    p1_z = start[2];
    p2_x = end[0];
    p2_y = end[1];
    p2_z = end[2];

    im_size_x = options->im_size_x;
    im_size_y = options->im_size_y;
    im_size_z = options->im_size_z;

	N_x=im_size_x+1;
	N_y=im_size_y+1;
	N_z=im_size_z+1;

	/* d: voxel size */
    d_x = options->d_x;
    d_y = options->d_y;
    d_z = options->d_z;

	/* b: grid offset from origin */
    b_x = options->b_x;
    b_y = options->b_y;
    b_z = options->b_z;



    /* use total lengh=alpha_max-alpha_min instead, to get everage, not sum. */
    /* moving back to original d_conv*/
    d_conv=sqrt( (p1_x-p2_x)*(p1_x-p2_x) + (p1_y-p2_y)*(p1_y-p2_y) + (p1_z-p2_z)*(p1_z-p2_z));


    x_defined =  !(equal_to_precision(p1_x,p2_x,PRECISION));
    y_defined =  !(equal_to_precision(p1_y,p2_y,PRECISION));
    z_defined =  !(equal_to_precision(p1_z,p2_z,PRECISION));
	
    if( !x_defined && !y_defined && !z_defined)
	return;

    if (x_defined) {
	alpha_x_min=min_dbl(alpha_fn(0, p1_x, p2_x, b_x, d_x), alpha_fn(N_x-1, p1_x, p2_x, b_x, d_x));
	alpha_x_max=max_dbl(alpha_fn(0, p1_x, p2_x, b_x, d_x), alpha_fn(N_x-1, p1_x, p2_x, b_x, d_x));
    }
    else {
	alpha_x_min=-2;
	alpha_x_max=2;
	i=(int) floor_j( phi(0.0, p1_x, p2_x, b_x, d_x));
	if ( i < 0 || i >= im_size_x)
	    return;
	alpha_x=2;
	i_min = 1;
	i_max = 0;
    }

    if(y_defined) {
	alpha_y_min=min_dbl(alpha_fn(0, p1_y, p2_y, b_y, d_y), alpha_fn(N_y-1, p1_y, p2_y, b_y, d_y));
	alpha_y_max=max_dbl(alpha_fn(0, p1_y, p2_y, b_y, d_y), alpha_fn(N_y-1, p1_y, p2_y, b_y, d_y));
    }
    else {
	alpha_y_min=-2;
	alpha_y_max=2;
	j=(int) floor_j( phi(0.0, p1_y, p2_y, b_y, d_y));
	if ( j < 0 || j >= im_size_y)
	    return;
	alpha_y=2;
	j_min = 1;
	j_max = 0;
    }

    		
    if(z_defined) {
	alpha_z_min=min_dbl(alpha_fn(0, p1_z, p2_z, b_z, d_z), alpha_fn(N_z-1, p1_z, p2_z, b_z, d_z));
	alpha_z_max=max_dbl(alpha_fn(0, p1_z, p2_z, b_z, d_z), alpha_fn(N_z-1, p1_z, p2_z, b_z, d_z));
    }
    else {
	alpha_z_min=-2;
	alpha_z_max=2;
	k=(int) floor_j( phi(0.0, p1_z, p2_z, b_z, d_z));
	if ( k < 0 || k >= im_size_z)
	    return;
	alpha_z=2;
	k_min = 1;
	k_max = 0;
    }
		
    alpha_min=max_dbl(0.0, max3_dbl(alpha_x_min, alpha_y_min, alpha_z_min));
    alpha_max=min_dbl(1.0, min3_dbl(alpha_x_max, alpha_y_max, alpha_z_max));

    /* if ray intersects voxel grid */
    if (alpha_min < alpha_max) {

	if (x_defined && p1_x < p2_x) {
	    if (equal_to_precision(alpha_min,alpha_x_min,PRECISION)==1)
		i_min=1;
	    else
		i_min = (int) ceil_j(phi(alpha_min, p1_x, p2_x, b_x, d_x));

	    if (equal_to_precision(alpha_max,alpha_x_max,PRECISION)==1)
		i_max = N_x - 1;
	    else
		i_max = (int) floor_j( phi(alpha_max, p1_x, p2_x, b_x, d_x));

	    alpha_x=alpha_fn(i_min, p1_x, p2_x, b_x, d_x);
	}

	else if (x_defined) {
	    if (equal_to_precision(alpha_min,alpha_x_min,PRECISION)==1)
		i_max=N_x-2;
	    else
		i_max = (int) floor_j(phi(alpha_min, p1_x, p2_x, b_x, d_x));

	    if (equal_to_precision(alpha_max,alpha_x_max,PRECISION)==1)
		i_min = 0;
	    else
		i_min = (int) ceil_j( phi(alpha_max, p1_x, p2_x, b_x, d_x));

	    alpha_x=alpha_fn(i_max, p1_x, p2_x, b_x, d_x);
	}


	if (y_defined && p1_y < p2_y) {
	    if (equal_to_precision(alpha_min,alpha_y_min,PRECISION)==1)
		j_min=1;
	    else
		j_min = (int) ceil_j(phi(alpha_min, p1_y, p2_y, b_y, d_y));


	    if (equal_to_precision(alpha_max, alpha_y_max,PRECISION)==1)
		j_max = N_y - 1;
	    else
		j_max = (int) floor_j( phi(alpha_max, p1_y, p2_y, b_y, d_y));

	    alpha_y=alpha_fn(j_min, p1_y, p2_y, b_y, d_y);
	}

	else if (y_defined) {

	    if (equal_to_precision(alpha_min,alpha_y_min,PRECISION)==1)
		j_max=N_y-2;
	    else
		j_max = (int) floor_j(phi(alpha_min, p1_y, p2_y, b_y, d_y));


	    if (equal_to_precision(alpha_max, alpha_y_max, PRECISION)==1)
		j_min = 0;
	    else
		j_min = (int) ceil_j( phi(alpha_max, p1_y, p2_y, b_y, d_y));

	    alpha_y=alpha_fn(j_max, p1_y, p2_y, b_y, d_y);
	}



	if (z_defined && p1_z < p2_z) {
	    if (equal_to_precision(alpha_min,alpha_z_min,PRECISION)==1)
		k_min=1;
	    else
		k_min = (int) ceil_j(phi(alpha_min, p1_z, p2_z, b_z, d_z));


	    if (equal_to_precision(alpha_max, alpha_z_max,PRECISION)==1)
		k_max = N_z - 1;
	    else
		k_max = (int) floor_j( phi(alpha_max, p1_z, p2_z, b_z, d_z));

	    alpha_z=alpha_fn(k_min, p1_z, p2_z, b_z, d_z);
	}

	else if (z_defined) {

	    if (equal_to_precision(alpha_min,alpha_z_min,PRECISION)==1)
		k_max=N_z-2;
	    else
		k_max = (int) floor_j(phi(alpha_min, p1_z, p2_z, b_z, d_z));


	    if (equal_to_precision(alpha_max, alpha_z_max, PRECISION)==1)
		k_min = 0;
	    else
		k_min = (int) ceil_j( phi(alpha_max, p1_z, p2_z, b_z, d_z));

	    alpha_z=alpha_fn(k_max, p1_z, p2_z, b_z, d_z);
	}


	N_p=(i_max - i_min +1) + (j_max - j_min + 1) + (k_max - k_min + 1);

	if (x_defined) {
	    i=(int) floor_j( phi( (min3_dbl(alpha_x, alpha_y, alpha_z) + alpha_min)/2, p1_x, p2_x, b_x, d_x) );
	alpha_x_u = d_x/fabs(p2_x-p1_x);
	}

	if (y_defined) {
	    j=(int) floor_j( phi( (min3_dbl(alpha_x, alpha_y, alpha_z) + alpha_min)/2, p1_y, p2_y, b_y, d_y) );
	alpha_y_u = d_y/fabs(p2_y-p1_y);
	}
	if (z_defined) {
	    k=(int) floor_j( phi( (min3_dbl(alpha_x, alpha_y, alpha_z) + alpha_min)/2, p1_z, p2_z, b_z, d_z) );
	alpha_z_u = d_z/fabs(p2_z-p1_z);
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

	for (n_count=1; n_count<N_p+1;n_count++) {


	    /* x smallest*/
	    if (x_defined && alpha_x <= alpha_y && alpha_x <= alpha_z) {
		/* ray intersects pixel(i,j) with length l_ij */

	      vol_data[ray_index] += (REAL)((alpha_x-alpha_c)*d_conv * ray_data);

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

	      vol_data[ray_index] += (REAL)((alpha_y-alpha_c)*d_conv * ray_data);

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

	      vol_data[ray_index] += (REAL)((alpha_z-alpha_c)*d_conv * ray_data);

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
	    
	    vol_data[ray_index] += (REAL)(l_ij * ray_data);
	}
	
    } /* of alpha_min < alpha_max */
    
    return;
}

static inline double r_i(const long i, const double r_0, const double step)
{
  // coord at grid index i
  return r_0 + i * step;
}

void backwardProjection(double *source_x, double *source_y, double *source_z,
			double *det_x, double *det_y, double *det_z,
			REAL *ray_data, REAL *vol_data, double *angles,
			struct jacobs_options *options, int n_angles,
			long n_rays_y, long n_rays_z, int size_z,
			double grid_offset[], double voxel_size[])
{
  int i, curr_angle, curr_ray_y, curr_ray_z;
  long ray_offset;
  double cos_curr_angle, sin_curr_angle;
  double start[3], end[3];
#pragma omp parallel shared(source_x, source_y, source_z, det_x, det_y, det_z, vol_data, angles, size_z, options, grid_offset, voxel_size) private(curr_angle, curr_ray_y, curr_ray_z, cos_curr_angle, sin_curr_angle, start, end, ray_offset), firstprivate(n_rays_y, n_rays_z, n_angles, ray_data)
  {
    struct jacobs_options local_opts = *options;
    int nz_offset = 0;
    int nz_step = 0;
    int nthreads = omp_get_num_threads();
    int threadid = omp_get_thread_num();
    int nzblocks = nthreads;
    int nz_size = size_z / nzblocks;
    int extra = size_z - nz_size * nzblocks;
    int half = (nthreads - 1) / 2;
    if (threadid <= half) {
      extra = (extra + 1) / 2;
      for (i = 0; i <= half; i++) {
	if (i > half - extra)
	  nz_step = nz_size + 1;
	else
	  nz_step = nz_size;
	if (i == threadid)
	  break;
	nz_offset += nz_step;
      }
    } else {
      nz_offset = size_z;
      extra = extra / 2;
      for (i = nthreads - 1; i > half; i--) {
	if (i <= half + extra)
	  nz_step = nz_size + 1;
	else
	  nz_step = nz_size;
	nz_offset -= nz_step;
	if (i == threadid)
	  break;
      }
    }

    if (nz_step > 0) {
      local_opts.im_size_z = nz_step;
      local_opts.b_z = grid_offset[2] + voxel_size[2] * nz_offset;

      for(curr_ray_z = 0; curr_ray_z < n_rays_z; curr_ray_z++) {
	double delta_z, alpha_z_0, alpha_z_N;
	double alpha_z_min, alpha_z_max, alpha_min, alpha_max;
	start[2] = *source_z;
	end[2] = det_z[curr_ray_z];

	delta_z = end[2] - start[2];
	alpha_z_0 = (r_i(0, local_opts.b_z, voxel_size[2]) - start[2]) / delta_z;
	alpha_z_N = (r_i(nz_step, local_opts.b_z, voxel_size[2]) - start[2]) / delta_z;
	alpha_z_min = min_dbl(alpha_z_0, alpha_z_N);
	alpha_z_max = max_dbl(alpha_z_0, alpha_z_N);
	alpha_min = max_dbl(0.0, alpha_z_min);
	alpha_max = min_dbl(1.0, alpha_z_max);

	if (alpha_min < alpha_max) {
	  for(curr_angle = 0; curr_angle < n_angles; curr_angle++) {
	    /* rotate source and detector positions by current angle */
	    cos_curr_angle = cos(angles[curr_angle]);
	    sin_curr_angle = sin(angles[curr_angle]);
        
	    start[0] = cos_curr_angle * (*source_x) - sin_curr_angle * (*source_y);
	    start[1] = sin_curr_angle * (*source_x) + cos_curr_angle * (*source_y);
        
	    ray_offset = curr_angle * n_rays_y * n_rays_z + curr_ray_z*n_rays_y;
        
	    /* loop over y values on detector */
	    for(curr_ray_y = 0; curr_ray_y < n_rays_y; curr_ray_y++) {
	      end[0] = cos_curr_angle * (*det_x) - sin_curr_angle * det_y[curr_ray_y];
	      end[1] = sin_curr_angle * (*det_x) + cos_curr_angle * det_y[curr_ray_y];
            
	      /* loop over z values on detector */
	      backproject_singledata(start, end,
				     (double)ray_data[ray_offset + curr_ray_y],
				     vol_data, &local_opts, nz_offset);
	    }
	  }
	}
      }
    }
  }    
}
