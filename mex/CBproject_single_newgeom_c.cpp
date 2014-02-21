/* function to perform cone beam back projection of Nikon XTek custom bay data */
/* uses an adaptation of David's Jacobs ray tracing code */

/* 06/09/2011 */

/* 07/09/2011 single precision version */

/* 11/11/2011 update for new geometry description */


#include <mex.h>
#include <cmath>
#include <matrix.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <omp.h>
#include "mex_types.hpp"
#include "instruments.hpp"
#include "project_line.hpp"
#include "cone_f.hpp"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int i, curr_angle, im_size, curr_ray_y, curr_ray_z;
	long n_rays_y, n_rays_z, n_angles, ray_offset;
	mwSize im_size_matlab[3];
	double *source_x, *source_y, *source_z, *det_x, *det_y, *det_z;
    double cos_curr_angle, sin_curr_angle;
	double start[3], end[3], *voxel_size, *grid_offset, *angles;
    float *ray_data, *vol_data;
    double *size_doubles;

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
    
    vol_data = (float *) mxGetData(prhs[9]);
    
    angles = mxGetPr(prhs[10]);
    
    n_angles = mxGetM(prhs[10]);
    
    if (mxGetM(prhs[9]) != (im_size_matlab[0] * im_size_matlab[1] * im_size_matlab[2]))
    {
        mexWarnMsgTxt("Mismatched number of voxels! Aborting.");
        return;
    }

  plhs[0] = mxCreateNumericMatrix(n_rays_y * n_rays_z * n_angles, 1, mxSINGLE_CLASS, mxREAL);
  ray_data = (float *) mxGetData(plhs[0]);

  CCPi::instrument::forward_project(*source_x, *source_y, *source_z, *det_x,
				    det_y, det_z, angles, ray_data,
				    (float *const) vol_data,
				    n_angles, n_rays_y, n_rays_z,
				    grid_offset, voxel_size,
				    im_size_matlab[0], im_size_matlab[1],
				    im_size_matlab[2]);
}
