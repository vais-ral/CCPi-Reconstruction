
#ifndef CCPI_PROJECT_LINE
#define CCPI_PROJECT_LINE

#include <cfloat>

namespace CCPi {

  template <class pixel_t, class voxel_t, bool backward>
  void generate_line(double L, double x, double y, double z,
		     const double L_x_inc, const double L_y_inc,
		     const double L_z_inc, const sl_int i_step,
		     const sl_int j_step, const sl_int k_step, sl_int ray_index,
		     pixel_t &pixel, voxel_t *const vol_data, const double tol1,
		     const double tol2);
  template <class pixel_t, class voxel_t, bool backward>
  void project_singledata(const real start[], const real end[],
			  pixel_t &ray_data, voxel_t *const vol_data,
			  const real b_x, const real b_y, const real b_z,
			  const real d_x, const real d_y, const real d_z,
			  const int im_size_x, const int im_size_y,
			  const int im_size_z, const sl_int z_offset);

}

#define PRECISION 0.00000001

template <class pixel_t, class voxel_t, bool backward>
void CCPi::generate_line(double L, double x, double y, double z,
			 const double L_x_inc, const double L_y_inc,
			 const double L_z_inc, const sl_int i_step,
			 const sl_int j_step, const sl_int k_step,
			 sl_int ray_index, pixel_t &pixel,
			 voxel_t *const vol_data, const double tol1,
			 const double tol2)
{
  double sum = 0.0;
  while (L > tol2) {
    if (x < y - tol1) {
      if (x < z - tol1) {
	L -= x;
	if (backward)
	  vol_data[ray_index] += voxel_t(x * pixel);
	else
	  sum += x * vol_data[ray_index];
	//i += i_inc;
	ray_index += i_step;
	y -= x;
	z -= x;
	x = L_x_inc;
      } else if (x > z - tol1 && x < z + tol1) {
	L -= x;
	if (backward)
	  vol_data[ray_index] += voxel_t(x * pixel);
	else
	  sum += x * vol_data[ray_index];
	//i += i_inc;
	//k += k_inc;
	ray_index += k_step + i_step;
	y -= x;
	z = L_z_inc;
	x = L_x_inc;
      } else {
	L -= z;
	if (backward)
	  vol_data[ray_index] += voxel_t(z * pixel);
	else
	  sum += z * vol_data[ray_index];
	//k += k_inc;
	ray_index += k_step;
	x -= z;
	y -= z;
	z = L_z_inc;
      }
    } else if (x > y - tol1 && x < y + tol1) {
      if (x < z) {
	L -= x;
	if (backward)
	  vol_data[ray_index] += voxel_t(x * pixel);
	else
	  sum += x * vol_data[ray_index];
	//i += i_inc;
	//j += j_inc;
	ray_index += j_step + i_step;
	y = L_y_inc;
	z -= x;
	x = L_x_inc;
      } else if (x > z - tol1 && x < z + tol1) {
	L -= x;
	if (backward)
	  vol_data[ray_index] += voxel_t(x * pixel);
	else
	  sum += x * vol_data[ray_index];
	//i += i_inc;
	//j += j_inc;
	//k += k_inc;
	ray_index += k_step + j_step + i_step;
	y = L_y_inc;
	z = L_z_inc;
	x = L_x_inc;
      } else {
	L -= z;
	if (backward)
	  vol_data[ray_index] += voxel_t(z * pixel);
	else
	  sum += z * vol_data[ray_index];
	//k += k_inc;
	ray_index += k_step;
	x -= z;
	y -= z;
	z = L_z_inc;
      }
    } else if (y < z - tol1) {
      L -= y;
      if (backward)
	vol_data[ray_index] += voxel_t(y * pixel);
      else
	sum += y * vol_data[ray_index];
      //j += j_inc;
      ray_index += j_step;
      x -= y;
      z -= y;
      y = L_y_inc;
    } else if (y > z - tol1 && y < z + tol1) {
      L -= y;
      if (backward)
	vol_data[ray_index] += voxel_t(y * pixel);
      else
	sum += y * vol_data[ray_index];
      //j += j_inc;
      //k += k_inc;
      ray_index += k_step + j_step;
      x -= y;
      z = L_z_inc;
      y = L_y_inc;
    } else {
      L -= z;
      if (backward)
	vol_data[ray_index] += voxel_t(z * pixel);
      else
	sum += z * vol_data[ray_index];
      //k += k_inc;
      ray_index += k_step;
      x -= z;
      y -= z;
      z = L_z_inc;
    }
  }
  if (!backward)
    pixel += pixel_t(sum);
}

template <class pixel_t, class voxel_t, bool backward>
void CCPi::project_singledata(const real start[], const real end[],
			      pixel_t &ray_data, voxel_t *const vol_data,
			      const real b_x, const real b_y, const real b_z,
			      const real d_x, const real d_y, const real d_z,
			      const int im_size_x, const int im_size_y,
			      const int im_size_z, const sl_int z_offset)
{
  const double source_x = start[0];
  const double source_y = start[1];
  const double source_z = start[2];
  const double detect_x = end[0];
  const double detect_y = end[1];
  const double detect_z = end[2];
  const sl_int nvoxels_x = im_size_x;
  const sl_int nvoxels_y = im_size_y;
  const sl_int nvoxels_z = im_size_z;
  const double voxel_size_x = d_x;
  const double voxel_size_y = d_y;
  const double voxel_size_z = d_z;
  const double voxel0_x = b_x;
  const double voxel0_y = b_y;
  const double voxel0_z = b_z;

  /* spacing between source and detector */
  const double delta_x = detect_x - source_x;
  const double delta_y = detect_y - source_y;
  const double delta_z = detect_z - source_z;

  const double delta_x_abs = std::abs(delta_x);
  const double delta_y_abs = std::abs(delta_y);
  const double delta_z_abs = std::abs(delta_z);

  if (delta_x_abs < PRECISION && delta_y_abs < PRECISION &&
      delta_z_abs < PRECISION)
    return;

  /* distance between source and detector */
  const double distance = sqrt(delta_x * delta_x + delta_y * delta_y
			     + delta_z * delta_z);
  /* coord at other edge of grid is */
  const double voxeln_x = voxel0_x + nvoxels_x * voxel_size_x;
  const double voxeln_y = voxel0_y + nvoxels_y * voxel_size_y;
  const double voxeln_z = voxel0_z + nvoxels_z * voxel_size_z;
  /* The line between the source and detector is r + alpha delta_r,
     where alpha = 0.0 at the source and 1.0 at the detector, so it
     intercepts the limits of the voxel grid at (r' - r) / delta_r */
  double alpha0_x;
  double alphan_x;
  if (delta_x_abs > PRECISION) {
    if (delta_x < 0.0) {
      alphan_x = (voxel0_x - source_x) / delta_x;
      alpha0_x = (voxeln_x - source_x) / delta_x;
    } else {
      alpha0_x = (voxel0_x - source_x) / delta_x;
      alphan_x = (voxeln_x - source_x) / delta_x;
    }
  } else {
    alpha0_x = -2.0;
    alphan_x = 2.0;
  }
  double alpha0_y;
  double alphan_y;
  if (delta_y_abs > PRECISION) {
    if (delta_y < 0.0) {
      alphan_y = (voxel0_y - source_y) / delta_y;
      alpha0_y = (voxeln_y - source_y) / delta_y;
    } else {
      alpha0_y = (voxel0_y - source_y) / delta_y;
      alphan_y = (voxeln_y - source_y) / delta_y;
    }
  } else {
    alpha0_y = -2.0;
    alphan_y = 2.0;
  }
  double alpha0_z;
  double alphan_z;
  if (delta_z_abs > PRECISION) {
    if (delta_z < 0.0) {
      alphan_z = (voxel0_z - source_z) / delta_z;
      alpha0_z = (voxeln_z - source_z) / delta_z;
    } else {
      alpha0_z = (voxel0_z - source_z) / delta_z;
      alphan_z = (voxeln_z - source_z) / delta_z;
    }
  } else {
    alpha0_z = -2.0;
    alphan_z = 2.0;
  }
  /* The line doesn't intercept the grid until all 3 coords intercept,
     so take max of alpha0, min of alphan */
  const double alpha0 = std::max(0.0, std::max(alpha0_x,
					       std::max(alpha0_y, alpha0_z)));
  const double alphan = std::min(1.0, std::min(alphan_x,
					       std::min(alphan_y, alphan_z)));
  /* if first intercept if beyond the second then it missed the grid */
  if (alpha0 < alphan) {
    const double dalpha = alphan - alpha0;
    /* x is the stride 1 index so order the loops so this is innermost,
       the alternative would be to order on the lengths in the different
       directions so the most steps in a direction is innermost which is
       a bit like Zhao and Reader (2003) */
    /* length of line within the grid */
    const double length = distance * dalpha;

    /* the delta alpha of a grid step is */
    const double dalpha_x = voxel_size_x / delta_x_abs;
    const double dalpha_y = voxel_size_y / delta_y_abs;
    const double dalpha_z = voxel_size_z / delta_z_abs;
    /* start and end coords in voxel grid */
    const double x_start = source_x + alpha0 * delta_x;
    const double y_start = source_y + alpha0 * delta_y;
    const double z_start = source_z + alpha0 * delta_z;
    /* converted into array indices */
    double x_start_idx = (x_start - voxel0_x) / voxel_size_x;
    double y_start_idx = (y_start - voxel0_y) / voxel_size_y;
    double z_start_idx = (z_start - voxel0_z) / voxel_size_z;

    // Length of line from one voxel to the next
    const double L_x_inc = dalpha_x * distance;
    const double L_y_inc = dalpha_y * distance;
    const double L_z_inc = dalpha_z * distance;

    const double tol1 = distance * 5.0 * DBL_EPSILON;

    double x_min;
    sl_int i_inc;
    double x;
    if (delta_x_abs < tol1) {
      x = distance;
      i_inc = 1;
      x_min = floor(x_start_idx + tol1);
    } else if (delta_x < 0.0) {
      x_min = floor(x_start_idx - tol1);
      i_inc = -1;
      if (x_start_idx > x_min - tol1 && x_start_idx < x_min + tol1)
	x = L_x_inc;
      else
	x = (x_start_idx - x_min) * L_x_inc;
    } else {
      x_min = floor(x_start_idx + tol1);
      i_inc = 1;
      if (x_start_idx > x_min - tol1 && x_start_idx < x_min + tol1)
	x = L_x_inc;
      else
	x = (1.0 - (x_start_idx - x_min)) * L_x_inc;
    }
    const sl_int i_min = (sl_int)x_min;
    sl_int i = i_min;
    double y_min;
    sl_int j_inc;
    double y;
    if (delta_y_abs < tol1) {
      y = distance;
      j_inc = 1;
      y_min = floor(y_start_idx + tol1);
    } else if (delta_y < 0.0) {
      y_min = floor(y_start_idx - tol1);
      j_inc = -1;
      if (y_start_idx > y_min - tol1 && y_start_idx < y_min + tol1)
	y = L_y_inc;
      else
	y = (y_start_idx - y_min) * L_y_inc;
    } else {
      y_min = floor(y_start_idx + tol1);
      j_inc = 1;
      if (y_start_idx > y_min - tol1 && y_start_idx < y_min + tol1)
	y = L_y_inc;
      else
	y = (1.0 - (y_start_idx - y_min)) * L_y_inc;
    }
    const sl_int j_min = (sl_int)y_min;
    sl_int j = j_min;
    double z_min;
    sl_int k_inc;
    double z;
    if (delta_z_abs < tol1) {
      z = distance;
      k_inc = 1;
      z_min = floor(z_start_idx + tol1);
    } else if (delta_z < 0.0) {
      z_min = floor(z_start_idx - tol1);
      k_inc = -1;
      if (z_start_idx > z_min - tol1 && z_start_idx < z_min + tol1)
	z = L_z_inc;
      else
	z = (z_start_idx - z_min) * L_z_inc;
    } else {
      z_min = floor(z_start_idx + tol1);
      k_inc = 1;
      if (z_start_idx > z_min - tol1 && z_start_idx < z_min + tol1)
	z = L_z_inc;
      else
	z = (1.0 - (z_start_idx - z_min)) * L_z_inc;
    }
    const sl_int k_min = (sl_int)z_min;
    sl_int k = k_min;

    double L = length;
    sl_int ray_index = (k+z_offset)*nvoxels_y*nvoxels_x + j*nvoxels_x + i;
    sl_int i_step = i_inc;
    sl_int j_step = j_inc * nvoxels_x;
    sl_int k_step = k_inc * nvoxels_y * nvoxels_x;
    // This needs to be big enough to allow for the accumulated rounding errors
    // of subtracting lengths from L, whilst not overrunning the edge of the
    // voxels. Otherwise we would also need to test i/j/k >=0 < nvoxels.
    const double tol2 = distance * 100.0 * DBL_EPSILON;
    // try to order so that most common test in loop is the first one.
    if (L_x_inc < L_y_inc + tol1) {
      if (L_x_inc < L_z_inc + tol1) {
	// x < y, z
	if (L_y_inc < L_z_inc + tol1)
	  // x, y, z
	  generate_line<pixel_t, voxel_t, backward>(L, x, y, z, L_x_inc, L_y_inc, L_z_inc, i_step, j_step,
			k_step, ray_index, ray_data, vol_data, tol1, tol2);
	else
	  // x, z, y
	  generate_line<pixel_t, voxel_t, backward>(L, x, z, y, L_x_inc, L_z_inc, L_y_inc, i_step, k_step,
			j_step, ray_index, ray_data, vol_data, tol1, tol2);
      } else
	generate_line<pixel_t, voxel_t, backward>(L, z, x, y, L_z_inc, L_x_inc, L_y_inc, k_step, i_step,
		      j_step, ray_index, ray_data, vol_data, tol1, tol2);
    } else {
      // x > y
      if (L_y_inc < L_z_inc + tol1) {
	if (L_x_inc < L_z_inc + tol1)
	  generate_line<pixel_t, voxel_t, backward>(L, y, x, z, L_y_inc, L_x_inc, L_z_inc, j_step, i_step,
			k_step, ray_index, ray_data, vol_data, tol1, tol2);
	else
	  generate_line<pixel_t, voxel_t, backward>(L, y, z, x, L_y_inc, L_z_inc, L_x_inc, j_step, k_step,
			i_step, ray_index, ray_data, vol_data, tol1, tol2);
      } else
	generate_line<pixel_t, voxel_t, backward>(L, z, y, x, L_z_inc, L_y_inc, L_x_inc, k_step, j_step,
		      i_step, ray_index, ray_data, vol_data, tol1, tol2);
    }
  }
}

#endif // CCPI_PROJECT_LINE
