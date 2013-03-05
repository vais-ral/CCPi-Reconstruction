
#ifndef CCPI_CONE_FORWARD
#define CCPI_CONE_FORWARD

template <class pixel_t, class voxel_t>
void CCPi::instrument::forward_project(const real source_x, const real source_y,
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
  int curr_angle, curr_ray_y, curr_ray_z;
  long ray_offset;
  real start[3], end[3];
#pragma omp parallel for shared(det_y, det_z, ray_data, phi) private(curr_angle, curr_ray_y, curr_ray_z, start, end, ray_offset), firstprivate(source_x, source_y, source_z, det_x, vol_data, n_angles, n_rays_y, n_rays_z, nx_voxels, ny_voxels, nz_voxels) schedule(dynamic)

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
      end[2] = sin_theta_angle * det_x + cos_theta_angle * det_z[curr_ray_z];
      xp = (cos_theta_angle * det_x - sin_theta_angle * det_z[curr_ray_z]);

      ray_offset = curr_angle * n_rays_y * n_rays_z + curr_ray_z*n_rays_y;

      /* loop over y values on detector */
      for(curr_ray_y = 0; curr_ray_y < n_rays_y; curr_ray_y++) {
	end[0] = cos_phi_angle * xp - sin_phi_angle * det_y[curr_ray_y];
	end[1] = sin_phi_angle * xp + cos_phi_angle * det_y[curr_ray_y];

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

#endif // CCPI_CONE_FORWARD
