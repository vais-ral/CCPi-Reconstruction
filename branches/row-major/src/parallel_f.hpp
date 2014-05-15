
#ifndef CCPI_PARALLEL_FORWARD
#define CCPI_PARALLEL_FORWARD

void CCPi::instrument::forward_project(const real_1d &det_y,
				       const real_1d &det_z,
				       const real_1d &phi, pixel_data &ray_data,
				       voxel_data &vol_data,
				       const int n_angles, const int n_rays_y,
				       const int n_rays_z,
				       const real grid_offset[3],
				       const real voxel_size[3],
				       const int nx_voxels,
				       const int ny_voxels,
				       const int nz_voxels)
{
  sl_int curr_angle, curr_ray_y, curr_ray_z;
  real cos_curr_angle, sin_curr_angle;
  real start[3], end[3];

  // set detector z to 2* the yz limits of the voxels, so it misses
  // longest voxel dim should be sqrt(3), so 2 should be safe
  real det_x = real(2.0) * std::max(std::abs(grid_offset[0]),
				    std::max(std::abs(grid_offset[1]),
					     std::abs(grid_offset[2])));

#pragma omp parallel for shared(det_y, det_z, ray_data, phi) private(curr_angle, curr_ray_y, curr_ray_z, start, end, cos_curr_angle, sin_curr_angle), firstprivate(det_x) schedule(dynamic)

  for(curr_ray_z = 0; curr_ray_z < n_rays_z; curr_ray_z++) {
    end[2] = det_z[curr_ray_z];
    start[2] = end[2];

    for(curr_angle = 0; curr_angle < n_angles; curr_angle++) {
      /* rotate source and detector positions by current angle */
      cos_curr_angle = std::cos(phi[curr_angle]);
      sin_curr_angle = std::sin(phi[curr_angle]);

      /* loop over y values on detector */
      for(curr_ray_y = 0; curr_ray_y < n_rays_y; curr_ray_y++) {
	end[0] = cos_curr_angle * det_x - sin_curr_angle * det_y[curr_ray_y];
	end[1] = sin_curr_angle * det_x + cos_curr_angle * det_y[curr_ray_y];
	start[0] = end[0] - real(3.0) * cos_curr_angle * det_x;
	start[1] = end[1] - real(3.0) * sin_curr_angle * det_x;

	/* loop over z values on detector */

	project_singledata<pixel_type, false>(start, end,
			   ray_data[curr_angle][curr_ray_y][curr_ray_z],
			   vol_data, grid_offset[0], grid_offset[1],
			   grid_offset[2], voxel_size[0], voxel_size[1],
			   voxel_size[2], nx_voxels, ny_voxels, nz_voxels,
						    nz_voxels, 0);
      }
    }
  }
}

#endif // CCPI_PARALLEL_FORWARD
