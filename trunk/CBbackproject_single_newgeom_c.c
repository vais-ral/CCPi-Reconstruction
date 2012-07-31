/* function to perform cone beam back projection of Nikon XTek custom bay data */
/* uses an adaptation of David's Jacobs ray tracing code */

/* 06/09/2011 */

/* 07/09/2011 single precision version */

/* 11/11/2011 update for new geometry description */


#include <mex.h>
#include <math.h>
#include <matrix.h>
#include <stdlib.h>
#include <string.h>
#include "jacobs_rays.h"
#include "omp.h"

extern
void backproject_singledata(const double start[], const double end[],
			    const double ray_data, float vol_data[],
			    const struct jacobs_options *options,
			    const long z_offset);

static inline double r_i(const long i, const double r_0, const double step)
{
  // coord at grid index i
  return r_0 + i * step;
}

static inline double min_dbl(double a, double b)
{
  return a < b ? a : b;
}

static inline double max_dbl(double a, double b)
{
  return a > b ? a : b;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	/* check number of arguments */
	if(nrhs != 11)
    {
		mexWarnMsgTxt("11 inputs required. Aborting");
		return;
    }
	else if(nlhs < 1)
    {
		mexWarnMsgTxt("At least one output required. Aborting");
		return;
	}



	int i, curr_angle, im_size, curr_ray_y, curr_ray_z, n_angles;
	long n_rays_y, n_rays_z, ray_offset;
	mwSize im_size_matlab[3];
	double *source_x, *source_y, *source_z, *det_x, *det_y, *det_z;
    double cos_curr_angle, sin_curr_angle;
	double start[3], end[3], *voxel_size, *grid_offset, *angles;
    float *ray_data, *vol_data;
    double *size_doubles;
    int nthreads, threadid;

	struct jacobs_options options;
	
	/* read variables */
	size_doubles = mxGetPr(prhs[0]);
	for(i = 0; i < 3; i++)
    {
	    im_size_matlab[i] = (int) size_doubles[i];
    }


    n_rays_y = mxGetM(prhs[5]);
    n_rays_z = mxGetM(prhs[6]);
    
 	source_x = mxGetPr(prhs[1]); 
 	source_y = mxGetPr(prhs[2]); 
 	source_z = mxGetPr(prhs[3]); 
 	det_x = mxGetPr(prhs[4]); 
 	det_y = mxGetPr(prhs[5]); 
 	det_z = mxGetPr(prhs[6]); 

	voxel_size = mxGetPr(prhs[7]);
	grid_offset = mxGetPr(prhs[8]);
    
    ray_data = (float *) mxGetData(prhs[9]);
    
    angles = mxGetPr(prhs[10]);
    
    n_angles = mxGetM(prhs[10]);
    
    if (mxGetM(prhs[9]) != (n_rays_y * n_rays_z * n_angles))
    {
        mexWarnMsgTxt("Mismatched number of data points! Aborting.");
        return;
    }

    plhs[0] = mxCreateNumericMatrix(im_size_matlab[0]*im_size_matlab[1]*im_size_matlab[2], 1, mxSINGLE_CLASS, mxREAL);
    vol_data = (float *) mxGetData(plhs[0]);

#pragma omp parallel shared(source_x, source_y, source_z, det_x, det_y, det_z, im_size, vol_data, angles, im_size_matlab, grid_offset, voxel_size) private(curr_angle, curr_ray_y, curr_ray_z, cos_curr_angle, sin_curr_angle, start, end, ray_offset, options, threadid, nthreads), firstprivate(n_rays_y, n_rays_z, n_angles, ray_data)
  {

    nthreads = omp_get_num_threads();
    threadid = omp_get_thread_num();

  int nzblocks = nthreads;
  int nz_size = im_size_matlab[2] / nzblocks;
  if (nz_size * nzblocks < im_size_matlab[2])
    nz_size++;
  int nz_offset = threadid * nz_size;
  int nz_step = nz_size;
    if (nz_offset + nz_step > im_size_matlab[2])
      nz_step = im_size_matlab[2] - nz_offset;
    printf("  z %d %d\n", nz_offset, nz_step);

    options.im_size_default = 0;
    options.im_size_x = im_size_matlab[0];
    options.im_size_y = im_size_matlab[1];
    options.im_size_z = nz_step;

    options.b_default = 0;
    options.b_x = grid_offset[0];
    options.b_y = grid_offset[1];
    options.b_z = grid_offset[2] + voxel_size[2] * nz_offset;

    options.d_x = voxel_size[0];
    options.d_y = voxel_size[1];
    options.d_z = voxel_size[2];

    for(curr_ray_z = 0; curr_ray_z < n_rays_z; curr_ray_z++) {
      start[2] = *source_z;
      end[2] = det_z[curr_ray_z];

      double delta_z = end[2] - start[2];
      double alpha_z_0 = (r_i(0, options.b_z, voxel_size[2]) - start[2]) / delta_z;
      double alpha_z_N = (r_i(nz_step, options.b_z, voxel_size[2]) - start[2]) / delta_z;
      double alpha_z_min = min_dbl(alpha_z_0, alpha_z_N);
      double alpha_z_max = max_dbl(alpha_z_0, alpha_z_N);
      double alpha_min = max_dbl(0.0, alpha_z_min);
      double alpha_max = min_dbl(1.0, alpha_z_max);

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
				   vol_data, &options, nz_offset);
	  }
	}
      }
    }
    //nz_offset += nz_step;
  }    
}










