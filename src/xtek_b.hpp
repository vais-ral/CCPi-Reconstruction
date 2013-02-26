
#ifndef CCPI_XTEK_BACKWARD
#define CCPI_XTEK_BACKWARD

static inline real r_i(const long i, const real r_0, const real step)
{
  // coord at grid index i
  return r_0 + i * step;
}

template <class pixel_t, class voxel_t>
void CCPi::instrument::backward_project(const real source_x,
					const real source_y,
					const real source_z, const real det_x,
					const real det_y[], const real det_z[],
					const real phi[], const real theta[],
					pixel_t ray_data[],
					voxel_t *const vol_data,
					const int n_angles, const int n_rays_y,
					const int n_rays_z,
					const real grid_offset[3],
					const real voxel_size[3],
					const int nx_voxels,
					const int ny_voxels,
					const int nz_voxels)
{
  int i, curr_angle, curr_ray_y, curr_ray_z;
  long ray_offset;
  real start[3], end[3];
#pragma omp parallel shared(det_y, det_z, phi, grid_offset, voxel_size) private(curr_angle, curr_ray_y, curr_ray_z, start, end, ray_offset), firstprivate(source_x, source_y, source_z, det_x, n_rays_y, n_rays_z, n_angles, ray_data, vol_data, nx_voxels, ny_voxels, nz_voxels)
  {
    int nz_offset = 0;
    int nz_step = 0;
    int nthreads = omp_get_num_threads();
    int threadid = omp_get_thread_num();
    int nzblocks = nthreads;
    int nz_size = nz_voxels / nzblocks;
    int extra = nz_voxels - nz_size * nzblocks;
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
      nz_offset = nz_voxels;
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
      //local_opts.im_size_z = nz_step;
      real b_z = grid_offset[2] + voxel_size[2] * nz_offset;

      for(curr_ray_z = 0; curr_ray_z < n_rays_z; curr_ray_z++) {

	for(curr_angle = 0; curr_angle < n_angles; curr_angle++) {
	  // rotate source and detector positions by current angle
	  // these are -ve of angles since they apply to the object and we
	  // are inverting the transformation onto the source/detector
	  // theta 0 is xy plane, +ve tilts object up wrt source -ve down
	  real cos_phi_angle = std::cos(phi[curr_angle]);
	  real sin_phi_angle = std::sin(phi[curr_angle]);
	  real cos_theta_angle = std::cos(theta[curr_angle]);
	  real sin_theta_angle = std::sin(theta[curr_angle]);
	  
	  real xp = cos_theta_angle * source_x - sin_theta_angle * source_z;
	  start[0] = cos_phi_angle * xp - sin_phi_angle * source_y;
	  start[1] = sin_phi_angle * xp + cos_phi_angle * source_y;
	  start[2] = sin_theta_angle * source_x + cos_theta_angle * source_z;
	  end[2] = sin_theta_angle * det_x
	    + cos_theta_angle * det_z[curr_ray_z];
	  xp = (cos_theta_angle * det_x - sin_theta_angle * det_z[curr_ray_z]);

	  real delta_z, alpha_z_0, alpha_z_N;
	  real alpha_z_min, alpha_z_max, alpha_min, alpha_max;

	  delta_z = end[2] - start[2];
	  alpha_z_0 = (r_i(0, b_z, voxel_size[2]) - start[2]) / delta_z;
	  alpha_z_N = (r_i(nz_step, b_z, voxel_size[2]) - start[2]) / delta_z;
	  alpha_z_min = std::min(alpha_z_0, alpha_z_N);
	  alpha_z_max = std::max(alpha_z_0, alpha_z_N);
	  alpha_min = std::max(0.0, alpha_z_min);
	  alpha_max = std::min(1.0, alpha_z_max);

	  if (alpha_min < alpha_max) {

	    ray_offset = curr_angle * n_rays_y * n_rays_z + curr_ray_z*n_rays_y;

	    /* loop over y values on detector */
	    for(curr_ray_y = 0; curr_ray_y < n_rays_y; curr_ray_y++) {
	      end[0] = cos_phi_angle * xp - sin_phi_angle * det_y[curr_ray_y];
	      end[1] = sin_phi_angle * xp + cos_phi_angle * det_y[curr_ray_y];

	      /* loop over z values on detector */
	      project_singledata<pixel_t, voxel_t, true>(start, end,
				 ray_data[ray_offset + curr_ray_y],
				 vol_data, grid_offset[0], grid_offset[1],
				 b_z, voxel_size[0], voxel_size[1],
				 voxel_size[2], nx_voxels, ny_voxels,
				 nz_step, nz_offset);
	    }
	  }
	}
      }
    }
  }
}

#endif // CCPI_XTEK_BACKWARD
