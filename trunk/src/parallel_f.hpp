
#ifndef CCPI_PARALLEL_FORWARD
#define CCPI_PARALLEL_FORWARD

template <class pixel_t, class voxel_t>
void CCPi::instrument::forward_project(const real det_y[], const real det_z[],
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
  int curr_angle, curr_ray_y, curr_ray_z;
  long ray_offset;
  real start[3], end[3];

  // set detector z to 2* the yz limits of the voxels, so it misses
  // longest voxel dim should be sqrt(3), so 2 should be safe
  real det_x = 2.0 * std::max(std::abs(grid_offset[0]),
			      std::max(std::abs(grid_offset[1]),
				       std::abs(grid_offset[2])));

#pragma omp parallel for shared(det_y, det_z, ray_data, phi) private(curr_angle, curr_ray_y, curr_ray_z, start, end, ray_offset), firstprivate(det_x) schedule(dynamic)

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

      end[2] = sin_theta_angle * det_x + cos_theta_angle * det_z[curr_ray_z];
      real xp = (cos_theta_angle * det_x - sin_theta_angle * det_z[curr_ray_z]);
      // dir is mapping of -1,0,0 direction from angles
      real dirz = - sin_theta_angle;
      start[2] = -end[2];
      real dirp = - cos_theta_angle;
      real dirx = cos_phi_angle * dirp;
      real diry = sin_phi_angle * dirp;

      // so the source is located in the direction of dir
      // at 2 * detx distance from the detector pixel
      // make it 3 * for safety
      start[2] = end[2] + 3.0 * det_x * dirz;
      ray_offset = curr_angle * n_rays_y * n_rays_z + curr_ray_z*n_rays_y;

      /* loop over y values on detector */
      for(curr_ray_y = 0; curr_ray_y < n_rays_y; curr_ray_y++) {
	end[0] = cos_phi_angle * xp - sin_phi_angle * det_y[curr_ray_y];
	end[1] = sin_phi_angle * xp + cos_phi_angle * det_y[curr_ray_y];
	start[0] = end[0] + 3.0 * det_x * dirx;
	start[1] = end[1] + 3.0 * det_x * diry;

	/* loop over z values on detector */

	project_singledata<pixel_t, voxel_t, false>(start, end,
			   ray_data[ray_offset + curr_ray_y],
			   vol_data, grid_offset[0], grid_offset[1],
			   grid_offset[2], voxel_size[0], voxel_size[1],
			   voxel_size[2], nx_voxels, ny_voxels, nz_voxels, 0);
      }
    }
  }
}

#endif // CCPI_PARALLEL_FORWARD
