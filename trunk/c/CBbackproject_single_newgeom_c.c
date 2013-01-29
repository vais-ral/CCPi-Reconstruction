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

extern void backwardProjection(double *source_x, double *source_y,
			       double *source_z, double *det_x, double *det_y,
			       double *det_z, float *ray_data, float *vol_data,
			       double *angles, struct jacobs_options *options,
			       int n_angles, long n_rays_y, long n_rays_z,
			       int size_z, double grid_offset[],
			       double voxel_size[]);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int i, curr_angle, im_size, curr_ray_y, curr_ray_z, n_angles;
	long n_rays_y, n_rays_z, ray_offset;
	mwSize im_size_matlab[3];
	double *source_x, *source_y, *source_z, *det_x, *det_y, *det_z;
    double cos_curr_angle, sin_curr_angle;
	double start[3], end[3], *voxel_size, *grid_offset, *angles;
    float *ray_data, *vol_data;
    double *size_doubles;

	struct jacobs_options options;
	
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

    options.im_size_default = 0;
    options.im_size_x = im_size_matlab[0];
    options.im_size_y = im_size_matlab[1];
    //options.im_size_z = nz_step;

    options.b_default = 0;
    options.b_x = grid_offset[0];
    options.b_y = grid_offset[1];
    //options.b_z = grid_offset[2] + voxel_size[2] * nz_offset;

    options.d_x = voxel_size[0];
    options.d_y = voxel_size[1];
    options.d_z = voxel_size[2];

    backwardProjection(source_x, source_y, source_z, det_x, det_y, det_z,
		       ray_data, vol_data, angles, &options, n_angles,
		       n_rays_y, n_rays_z, im_size_matlab[2], grid_offset,
		       voxel_size);
}
