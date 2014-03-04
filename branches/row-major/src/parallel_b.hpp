
#ifndef CCPI_PARALLEL_BACKWARD
#define CCPI_PARALLEL_BACKWARD

static inline real pb_r_i(const sl_int i, const real r_0, const real step)
{
  // coord at grid index i
  return r_0 + i * step;
}

template <class pixel_t, class voxel_t>
void CCPi::instrument::backward_project(const real det_y[], const real det_z[],
					const real phi[], pixel_t ray_data[],
					voxel_t *const vol_data,
					const int n_angles, const int n_rays_y,
					const int n_rays_z,
					const real grid_offset[3],
					const real voxel_size[3],
					const int nx_voxels,
					const int ny_voxels,
					const int nz_voxels)
{
  sl_int curr_angle, curr_ray_y, curr_ray_z;
  sl_int ray_offset;
  real cos_curr_angle, sin_curr_angle;
  real start[3], end[3];

  real det_x = real(2.0) * std::max(std::abs(grid_offset[0]),
				    std::max(std::abs(grid_offset[1]),
					     std::abs(grid_offset[2])));

#pragma omp parallel shared(det_y, det_z, phi, grid_offset, voxel_size) private(curr_angle, curr_ray_y, curr_ray_z, start, end, ray_offset, cos_curr_angle, sin_curr_angle), firstprivate(det_x)
  {
    int nz_offset = 0;
    int nz_step = 0;
    int nthreads = omp_get_num_threads();
    int threadid = omp_get_thread_num();
    int nzblocks = nthreads;
    int nz_size = nz_voxels / nzblocks;
    int extra = nz_voxels - nz_size * nzblocks;
    if (threadid < extra) {
      nz_step = nz_size + 1;
      nz_offset = threadid * nz_step;
    } else {
      nz_step = nz_size;
      nz_offset = threadid * nz_size + extra;
    }

    if (nz_step > 0) {
      //local_opts.im_size_z = nz_step;
      real b_z = grid_offset[2] + voxel_size[2] * nz_offset;
      // start == end so delta_z = 0.0
      real min_z = pb_r_i(0, b_z, voxel_size[2]);
      real max_z = pb_r_i(nz_step, b_z, voxel_size[2]);

      for (curr_ray_z = 0; curr_ray_z < n_rays_z; curr_ray_z++) {
	end[2] = det_z[curr_ray_z];
	start[2] = end[2];

	// make sure the line is with the voxel box range
	if (min_z <= end[2] and end[2] < max_z) {
	  for(curr_angle = 0; curr_angle < n_angles; curr_angle++) {
	    /* rotate source and detector positions by current angle */
	    cos_curr_angle = std::cos(phi[curr_angle]);
	    sin_curr_angle = std::sin(phi[curr_angle]);

	    ray_offset = curr_angle * sl_int(n_rays_y) * sl_int(n_rays_z)
	      + curr_ray_z;

	    /* loop over y values on detector */
	    for(curr_ray_y = 0; curr_ray_y < n_rays_y; curr_ray_y++) {
	      end[0] = cos_curr_angle * det_x
		- sin_curr_angle * det_y[curr_ray_y];
	      end[1] = sin_curr_angle * det_x
		+ cos_curr_angle * det_y[curr_ray_y];
	      start[0] = end[0] - real(3.0) * cos_curr_angle * det_x;
	      start[1] = end[1] - real(3.0) * sin_curr_angle * det_x;

	      /* loop over z values on detector */
	      project_singledata<pixel_t, voxel_t, true>(start, end,
			   ray_data[ray_offset + curr_ray_y * sl_int(n_rays_z)],
				 vol_data, grid_offset[0], grid_offset[1],
				 b_z, voxel_size[0], voxel_size[1],
				 voxel_size[2], nx_voxels, ny_voxels,
				 nz_voxels, nz_step, nz_offset);
	    }
	  }
	}
      }
    }
  }
}

#endif // CCPI_PARALLEL_BACKWARD
