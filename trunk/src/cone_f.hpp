
#ifndef CCPI_CONE_FORWARD
#define CCPI_CONE_FORWARD

template <class pixel_t>
void CCPi::instrument::forward_project(const real source_x, const real source_y,
				       const real source_z, const real det_x,
				       const real_1d &det_y,
				       const real_1d &det_z,
				       const real_1d &phi, pixel_t ray_data[],
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
  sl_int ray_offset;
  real cos_curr_angle, sin_curr_angle;
  real start[3], end[3];
#pragma omp parallel for shared(det_y, det_z, ray_data, phi) private(curr_angle, curr_ray_y, curr_ray_z, start, end, ray_offset, cos_curr_angle, sin_curr_angle) schedule(dynamic)

  for(curr_ray_z = 0; curr_ray_z < n_rays_z; curr_ray_z++) {
    start[2] = source_z;
    end[2] = det_z[curr_ray_z];

    for(curr_angle = 0; curr_angle < n_angles; curr_angle++) {
      /* rotate source and detector positions by current angle */
      cos_curr_angle = std::cos(phi[curr_angle]);
      sin_curr_angle = std::sin(phi[curr_angle]);

      start[0] = cos_curr_angle * source_x - sin_curr_angle * source_y;
      start[1] = sin_curr_angle * source_x + cos_curr_angle * source_y;

      ray_offset = curr_angle * sl_int(n_rays_y) * sl_int(n_rays_z)
	+ curr_ray_z * sl_int(n_rays_y);

      /* loop over y values on detector */
      for(curr_ray_y = 0; curr_ray_y < n_rays_y; curr_ray_y++) {
	end[0] = cos_curr_angle * det_x - sin_curr_angle * det_y[curr_ray_y];
	end[1] = sin_curr_angle * det_x + cos_curr_angle * det_y[curr_ray_y];

	/* loop over z values on detector */

	project_singledata<pixel_t, false>(start, end,
			   ray_data[ray_offset + curr_ray_y],
			   vol_data, grid_offset[0], grid_offset[1],
			   grid_offset[2], voxel_size[0], voxel_size[1],
			   voxel_size[2], nx_voxels, ny_voxels, nz_voxels, 0);
      }
    }
  }
}

#endif // CCPI_CONE_FORWARD
