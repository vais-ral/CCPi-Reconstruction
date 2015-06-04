
#include <iostream>
#include <cmath>
#include "base_types.hpp"
#include "filters.hpp"
#include "instruments.hpp"
#ifdef HAS_NEXUS
#  include "nexus.hpp"
#endif // HAS_NEXUS
#include "utils.hpp"
#include "ui_calls.hpp"

bool CCPi::Diamond::setup_experimental_geometry(const std::string path,
						const std::string file,
						const real rotation_centre,
						const int pixels_per_voxel,
						const bool phantom)
{
  if (phantom)
    return create_phantom();
  else {
    name = file;
    return read_data_size(path, rotation_centre, pixels_per_voxel);
  }
}

bool CCPi::Diamond::setup_experimental_geometry(const numpy_3d &pix_array,
						const numpy_1d &angle_array,
						const real rotation_centre,
						const int pixels_per_voxel)
{
  bool ok = true;
  int nangles = (int)angle_array.shape()[0];
  if (nangles < 1) {
    report_error("Bad angle array");
    ok = false;
  } else {
    const pixel_data::size_type *s = pix_array.shape();
    int na = (int)s[0];
    int nv_pixels = (int)s[1];
    int nh_pixels = (int)s[2];
    if (na != nangles) {
      report_error("Number of projections doesn't match angle array");
      ok = false;
    } else if (nh_pixels < 1 or nv_pixels < 1) {
      report_error("Bad array index for pixels");
      ok = false;
    } else {
      // dummy
      real hsize = 1.0;
      real vsize = 1.0;
      // copied from read_data_size
      nv_pixels = calc_v_alignment(nv_pixels, pixels_per_voxel, false);
      real shift = 0.0;
      if (rotation_centre > 0.0) {
	// shift h_pixels
	// if 4 pixels and centre = 2.0 nothing happens
	real pixel_shift = rotation_centre - real(nh_pixels) / 2.0;
	shift = hsize * pixel_shift;
      }
      // store data in class
      real_1d &h_pixels = set_h_pixels(nh_pixels);
      h_pixels[0] = - ((nh_pixels - 1) * hsize) / real(2.0) - shift;
      for (int i = 1; i < nh_pixels; i++)
	h_pixels[i] = h_pixels[0] + real(i) * hsize;
      real_1d &v_pixels = set_v_pixels(nv_pixels);
      v_pixels[0] = - ((nv_pixels - 1) * vsize) / real(2.0);
      for (int i = 1; i < nv_pixels; i++)
	v_pixels[i] = v_pixels[0] + real(i) * vsize;
      real_1d &angles = set_phi(nangles);
      for (int i = 0; i < nangles; i++)
	angles[i] = M_PI * (angle_array[i] / 180.0);
    }
  }
  return ok;
}

bool CCPi::Diamond::setup_experimental_geometry(const numpy_3d &pix_array,
						const numpy_1d &angle_array,
						const numpy_1d &h_offsets,
						const numpy_1d &v_offsets,
						const int pixels_per_voxel,
						const real source_x,
						const real detector_x,
						const real pixel_h_size,
						const real pixel_v_size)
{
  report_error("Diamond Avizo interface not implemented");
  return false;
}

bool CCPi::Diamond::create_phantom()
{
  //set_source(-250.0, 0.0, 0.0);
  //set_detector(737.0);
  int c;
  for (c = 0; c < 1000; c++) {
    //real p = -100.0 + c * 0.250;
    real p = -99.8046875 + real(c) * 0.390625;
    if (p >= 100.0 + 0.001)
      break;
  }
  real_1d &h_pixels = set_h_pixels(c);
  for (int i = 0; i < c; i++) {
    //real p = -100.0 + i * 0.250;
    h_pixels[i] = -99.8046875 + real(i) * 0.390625;
  }
  for (c = 0; c < 1000; c++) {
    //real p = -100.0 + real(c) * 0.250;
    real p = -99.8046875 + real(c) * 0.390625;
    if (p >= 100.0 + 0.001)
      break;
  }
  real_1d &v_pixels = set_v_pixels(c);
  for (int i = 0; i < c; i++) {
    //real p = -100.0 + real(i) * 0.250;
    v_pixels[i] = -99.8046875 + real(i) * 0.390625;
  }
  // 501 values from 0 to 2pi
  //geom.angles = linspace(0,2*pi,501);
  // lose the 2pi which is a duplicate of 0
  //geom.angles = geom.angles(1:500);
  const int nangles = 1000;
  real_1d &pangles = set_phi(nangles);
  real step = M_PI / real(nangles);
  for (int i = 0; i < nangles; i++)
    pangles[i] = real(i) * step;
  return true;
}

bool CCPi::Diamond::read_data_size(const std::string path,
				   const real rotation_centre,
				   const int pixels_per_voxel)
{
  bool ok = false;
  std::string fullname;
  combine_path_and_name(path, name,fullname);
  pixel_data pixels(boost::extents[1][1][1]);
  int nh_pixels = 0;
  int nv_pixels = 0;
  std::vector<real> angles;
  int nangles = 0;
  real hsize = 0.0;
  real vsize = 0.0;
  pixel_2d i_dark(boost::extents[1][1]);
  pixel_2d f_dark(boost::extents[1][1]);
  pixel_2d i_bright(boost::extents[1][1]);
  pixel_2d f_bright(boost::extents[1][1]);
#ifdef HAS_NEXUS
  ok = read_NeXus(pixels, i_dark, f_dark, i_bright, f_bright, nh_pixels,
		  nv_pixels, angles, nangles, hsize, vsize, fullname,
		  false, false, 0, 10000, 0);
#endif // HAS_NEXUS
  if (ok) {
    nv_pixels = calc_v_alignment(nv_pixels, pixels_per_voxel, false);
    real shift = 0.0;
    if (rotation_centre > 0.0) {
      // shift h_pixels
      // if 4 pixels and centre = 2.0 nothing happens
      real pixel_shift = rotation_centre - real(nh_pixels) / 2.0;
      shift = hsize * pixel_shift;
    }
    // store data in class
    real_1d &h_pixels = set_h_pixels(nh_pixels);
    h_pixels[0] = - ((nh_pixels - 1) * hsize) / real(2.0) - shift;
    for (int i = 1; i < nh_pixels; i++)
      h_pixels[i] = h_pixels[0] + real(i) * hsize;
    real_1d &v_pixels = set_v_pixels(nv_pixels);
    v_pixels[0] = - ((nv_pixels - 1) * vsize) / real(2.0);
    for (int i = 1; i < nv_pixels; i++)
      v_pixels[i] = v_pixels[0] + real(i) * vsize;
    real_1d &phi = set_phi(nangles);
    for (int i = 0; i < nangles; i++)
      phi[i] = angles[i];
  }
  return ok;
}

bool CCPi::Diamond::read_scans(const std::string path, const int offset,
			       const int block_size, const bool first,
			       const bool phantom)
{
  if (phantom)
    return build_phantom(offset, block_size);
  else
    return read_data(path, offset, block_size, first);
}

bool CCPi::Diamond::build_phantom(const int offset, const int block_size)
{
  // build small voxel set to project to produce initial pixels
  int nx = 512;
  int ny = 512;
  int nz = 512;

  const real_1d &h_pixels = get_h_pixels();
  const real_1d &v_pixels = get_all_v_pixels();
  int nh = get_num_h_pixels();
  int nv = total_num_v_pixels();
  real hmin = std::abs(h_pixels[0]);
  real hmax = h_pixels[nh - 1] + (h_pixels[1] - h_pixels[0]);
  real vmin = std::abs(v_pixels[0]);
  real vmax = v_pixels[nv - 1] + (v_pixels[1] - v_pixels[0]);
  real hlim = std::max(hmax, hmin);
  real vlim = std::max(vmax, vmin);
  real lim = std::max(hlim, vlim);
  real voxel_size[3];
  voxel_size[0] = real(2.0) * lim / real(nx);
  voxel_size[1] = voxel_size[0];
  voxel_size[2] = voxel_size[0];
  real image_vol[3];
  image_vol[0] = voxel_size[0] * nx;
  image_vol[1] = voxel_size[1] * ny;
  image_vol[2] = voxel_size[2] * nz;
  real image_offset[3];
  image_offset[0] = -image_vol[0] / 2;
  image_offset[1] = -image_vol[1] / 2;
  image_offset[2] = -image_vol[2] / 2;

  // set up phantom volume
  voxel_data x(boost::extents[nx][ny][nz], boost::c_storage_order());
  //sl_int n_vox = nx * ny * nz;
  for (sl_int k = 0; k < nx; k++)
    for (sl_int j = 0; j < ny; j++)
      for (sl_int i = 0; i < nz; i++)
	x[k][j][i] = 0.0;

  // add cubes - column major
  for (int i = 108-1; i < 189; i++) {
    for (int j = 108-1; j < 189; j++) {
      for (int k = 58-1; k < 139; k++) {
        x[i][j][k] = 1;
      }
    }
  }
  for (int i = 190-1; i < 271; i++) {
    for (int j = 190-1; j < 271; j++) {
      for (int k = 140-1; k < 221; k++) {
        x[i][j][k] = 1;
      }
    }
  }
  for (int i = 272-1; i < 353; i++) {
    for (int j = 272-1; j < 353; j++) {
      for (int k = 222-1; k < 303; k++) {
        x[i][j][k] = 1;
      }
    }
  }

  set_v_offset(offset);
  pixel_data &pixels = create_pixel_data();
  // perform projection step
  safe_forward_project(pixels, x, image_offset, voxel_size, nx, ny, nz);
  //delete [] x;
  return true;
}

bool CCPi::Diamond::read_data(const std::string path, const int offset,
			      const int block_size, const bool first)
{
  bool ok = false;
  std::string fullname;
  combine_path_and_name(path, name,fullname);
#ifdef HAS_NEXUS
  real hsize = 0.0;
  real vsize = 0.0;
  int nh_pixels = 0;
  int nv_pixels = 0;
#endif // HAS_NEXUS
  std::vector<real> angles(2 * get_num_angles());
  int nangles = 0;
  int nv = get_num_v_pixels();
  int nh = get_num_h_pixels();
  if (first)
    (void) create_pixel_data();
  pixel_data &pixels = get_pixel_data();
  pixel_2d i_dark(boost::extents[nh][nv]);
  pixel_2d f_dark(boost::extents[nh][nv]);
  pixel_2d i_bright(boost::extents[nh][nv]);
  pixel_2d f_bright(boost::extents[nh][nv]);
  for (sl_int i = 0; i < nh; i++) {
    for (sl_int j = 0; j < nv; j++) {
      i_dark[i][j] = 0.0;
      f_dark[i][j] = 0.0;
      i_bright[i][j] = 0.0;
      f_bright[i][j] = 0.0;
    }
  }
  int v_offset = get_data_v_offset();
  int v_size = get_data_v_size();
  int v_end = v_offset + v_size;
#ifdef HAS_NEXUS
  int bsize = block_size;
  // need to be careful that we only read bits for the data set range
  if (block_size > v_size)
    bsize = v_size;
  else if (offset + bsize > v_end)
    bsize = v_end - offset;
  else if (offset == 0)
    bsize = block_size - v_offset;
  if (bsize > 0)
    ok = read_NeXus(pixels, i_dark, f_dark, i_bright, f_bright, nh_pixels,
		    nv_pixels, angles, nangles, hsize, vsize, fullname,
		    false, true, offset, bsize, v_offset);
  else // Its an empty data range
    ok = true;
#endif // HAS_NEXUS
  nangles = get_num_angles();
  // store data in class - now done by read_data_sizes
  set_v_offset(offset);
  if (ok) {
    const real_1d &angles = get_phi();
    // linear interpolate bright/dark frames. Todo - something else?
    // Todo, also use initial final bright/dark angles? rather than assuming
    // initial/final sample angles?
    pixel_2d dark(boost::extents[nh][nv]);
    pixel_2d bright(boost::extents[nh][nv]);
    for (int i = 0; i < nangles; i++) {
      // Based on fbp code, interpolate bright/dark
      real w = (angles[i] - angles[0]) / (angles[nangles - 1] - angles[0]);
      for (sl_int k = 0; k < nh; k++)
	for (sl_int j = v_offset; j < v_end; j++)
	  dark[k][j] = i_dark[k][j] * (real(1.0) - w) + f_dark[k][j] * w;
      for (sl_int k = 0; k < nh; k++)
	for (sl_int j = v_offset; j < v_end; j++)
	  bright[k][j] = i_bright[k][j] * (real(1.0) - w) + f_bright[k][j] * w;
      // subtract dark from data/bright
      // and clamp min data/bright value to 0.1
      for (sl_int k = 0; k < nh; k++) {
	for (sl_int j = v_offset; j < v_end; j++) {
	  bright[k][j] -= dark[k][j];
	  if (bright[k][j] < real(0.1))
	    bright[k][j] = 0.1;
	}
      }
      for (sl_int k = 0; k < nh; k++) {
	for (sl_int j = v_offset; j < v_end; j++) {
	  pixels[i][k][j] -= dark[k][j];
	  if (pixels[i][k][j] < real(0.1))
	    pixels[i][k][j] = 0.1;
	}
      }
      // scale each data pixel by bright pixel
      for (sl_int k = 0; k < nh; k++)
	for (sl_int j = v_offset; j < v_end; j++)
	  pixels[i][k][j] /= bright[k][j];
      // clamp to 1.0
      for (sl_int k = 0; k < nh; k++) {
	for (sl_int j = v_offset; j < v_end; j++) {
	  if (pixels[i][k][j] > real(1.0))
	    pixels[i][k][j] = 1.0;
	}
      }
    }
    // take -ve log, due to exponential extinction in sample.
    for (int i = 0; i < nangles; i++) {
      for (sl_int k = 0; k < nh; k++) {
	// initialise aligned bits that don't contain data
	for (sl_int j = 0; j < v_offset; j++) {
	  pixels[i][k][j] = 0.0;
	}
	for (sl_int j = v_offset; j < v_end; j++)
	  pixels[i][k][j] = - std::log(pixels[i][k][j]);
	for (sl_int j = v_end; j < nv; j++) {
	  pixels[i][k][j] = 0.0;
	}
      }
    }
    //find_centre(get_num_v_pixels() / 2 + 1);
    if (use_high_peaks)
      high_peaks_before(hp_jump, hp_num_pix);
    if (use_ring_artefacts)
      ring_artefact_removal(ra_algorithm, aml_param_n, aml_param_r,
			    ra_num_series);
    //if (use_intensity_norm) - port fbp intensity_norm routine
    //normalise_intensity();
  }
  return ok;
}

bool CCPi::Diamond::read_scans(const numpy_3d &pixel_array,
			       const int offset, const int block_size)
{
  // bits from read_data
  int nv = get_num_v_pixels();
  int nh = get_num_h_pixels();
  int nangles = get_num_angles();
  int v_offset = get_data_v_offset();
  int v_size = get_data_v_size();
  int v_end = v_offset + v_size;
  /*
  int bsize = block_size;
  if (block_size > v_size)
    bsize = v_size;
  else if (offset + bsize > v_end)
    bsize = v_end - offset;
  else if (offset == 0)
    bsize = block_size - v_offset;
  */
  set_v_offset(offset);
  // Copy pixels
  pixel_data &pixels = create_pixel_data();
  for (int i = 0; i < nangles; i++) {
    for (sl_int j = 0; j < nh; j++) {
      // initialise aligned bits that don't contain data
      for (sl_int k = 0; k < v_offset; k++) {
	pixels[i][j][k] = 0.0;
      }
      for (sl_int k = v_offset; k < v_end; k++)
	pixels[i][j][k] = - std::log(pixel_array[i][k - v_offset][j]);
      for (sl_int k = v_end; k < nv; k++) {
	pixels[i][j][k] = 0.0;
      }
    }
  }
  return true;
}

bool CCPi::Diamond::finish_voxel_geometry(real voxel_origin[3],
					  real voxel_size[3], const int nx,
					  const int ny, const int nz) const
{
  /*
  const real_1d &h_pixels = get_h_pixels();
  const real_1d &v_pixels = get_all_v_pixels();
  int nh = get_num_h_pixels();
  int nv = total_num_v_pixels();
  real hrange = std::abs(h_pixels[nh - 1] - h_pixels[0]);
  real hsize = hrange / (nh - 1);
  real vrange = std::abs(v_pixels[nv - 1] - v_pixels[0]);
  real vsize = vrange / (nv - 1);
  int hs = std::min(nx, ny);
  hsize *= (real(nh) / real(hs));
  vsize *= (real(nv) / nz);
  */
  voxel_size[0] = h_vox_size;
  voxel_size[1] = h_vox_size;
  voxel_size[2] = v_vox_size;
  voxel_origin[0] = -voxel_size[0] * real(nx) / real(2.0); // + offset[0];
  voxel_origin[1] = -voxel_size[1] * real(ny) / real(2.0); // + offset[1];
  voxel_origin[2] = -voxel_size[2] * real(nz) / real(2.0); // + offset[2];
  return true;
}

void CCPi::Diamond::get_xy_size(int &nx, int &ny, const int pixels_per_voxel)
{
  const real_1d &h_pixels = get_h_pixels();
  int nh = get_num_h_pixels();
  real hrange = h_pixels[1] - h_pixels[0];
  real hmax = std::max(std::abs(h_pixels[0]), std::abs(h_pixels[nh - 1]));
  int n = int(std::ceil((2.0 * hmax) / hrange));
  nx = n / pixels_per_voxel;
  if (n % pixels_per_voxel != 0)
    nx++;
  h_vox_size = hrange * real(pixels_per_voxel);
  // check that we've covered a big enough region
  real xmin = -h_vox_size * real(nx) / real(2.0);
  real xmax =  h_vox_size * real(nx) / real(2.0);
  while (xmin >= h_pixels[0] or xmax <= h_pixels[nh - 1]) {
    nx++;
    xmin = -h_vox_size * real(nx) / real(2.0);
    xmax =  h_vox_size * real(nx) / real(2.0);
  }
  ny = nx;
  const real_1d &v_pixels = get_all_v_pixels();
  v_vox_size = (v_pixels[1] - v_pixels[0]) * real(pixels_per_voxel);
}

void CCPi::Diamond::apply_beam_hardening()
{
  // Todo - does this belong in the base class?
  pixel_data &pixels = get_pixel_data();
  for (sl_int i = 0; i < get_num_angles(); i++)
    for (sl_int k = 0; k < get_num_h_pixels(); k++)
      for (sl_int j = 0; j < get_num_v_pixels(); j++)
	pixels[i][k][j] = pixels[i][k][j] * pixels[i][k][j];
}

void CCPi::Diamond::high_peaks_before(const real jump, const int num_pix)
{
  // ported from high_peaks_before routine in Manchester/Diamond FBP code
  sl_int nh = get_num_h_pixels();
  sl_int nv = get_num_v_pixels();
  sl_int nangles = get_num_angles();
  if (jump < 0.0) {
    report_error("High peaks - negative jump");
    return;
  }
  if (num_pix < 1 or num_pix >= nangles) {
    report_error("High peaks - wrong number of neighbours");
    return;
  }
  pixel_data &pixels = get_pixel_data();
  // correct each slice
  pixel_2d temp1(boost::extents[nangles - 1][nh]);
  //pixel_2d temp2(boost::extents[nangles - 1][nh]);
  // Todo - better loop structure - might require larger temp space.
  for (int v = 0; v < nv; v++) {
    // iterate num_pix times
    for (int k = 0; k < num_pix; k++) {
      // ToptoBottom
      //ippsSub
      for (int a = 0; a < nangles - 1; a++) {
	for (int h = 0; h < nh; h++)
	  temp1[a][h] = pixels[a][h][v] - pixels[a + 1][h][v];
      }
      /*
      // ippsThreshold_LTValGTVal
      for (int a = 0; a < nangles - 1; a++) {
	for (int h = 0; h < nh; h++) {
	  if (temp1[a][h] < -jump)
	    temp2[a][h] = 1.0;
	  else if (temp2[a][h] > -jump)
	    temp2[a][h] = 0.0;
	  else
	    temp2[a][h] = temp1[a][h];
	}
      }
      */
      // ippsTheshold_LT
      for (int a = 0; a < nangles - 1; a++) {
	for (int h = 0; h < nh; h++) {
	  if (temp1[a][h] < -jump)
	    temp1[a][h] = -jump;
	}
      }
      /*
      // ippsSum - could probably do this with LTValGTVal
      pixel_type sum = 0.0;
      for (int a = 0; a < nangles - 1; a++) {
	for (int h = 0; h < nh; h++) {
	  sum += temp2[a][h];
	}
      }
      */
      // ippsAdd
      for (int a = 0; a < nangles - 1; a++) {
	for (int h = 0; h < nh; h++)
	  pixels[a][h][v] = pixels[a + 1][h][v] + temp1[a][h];
      }
      //BottomtoTop
      //ippsSub
      for (int a = 0; a < nangles - 1; a++) {
	for (int h = 0; h < nh; h++)
	  temp1[a][h] = pixels[a + 1][h][v] - pixels[a][h][v];
      }
      /*
      // ippsThreshold_LTValGTVal
      for (int a = 0; a < nangles - 1; a++) {
	for (int h = 0; h < nh; h++) {
	  if (temp1[a][h] < -jump)
	    temp2[a][h] = 1.0;
	  else if (temp2[a][h] > -jump)
	    temp2[a][h] = 0.0;
	  else
	    temp2[a][h] = temp1[a][h];
	}
      }
      */
      // ippsTheshold_LT
      for (int a = 0; a < nangles - 1; a++) {
	for (int h = 0; h < nh; h++) {
	  if (temp1[a][h] < -jump)
	    temp1[a][h] = -jump;
	}
      }
      /*
      // ippsSum - could probably do this with LTValGTVal
      pixel_type sum = 0.0;
      for (int a = 0; a < nangles - 1; a++) {
	for (int h = 0; h < nh; h++) {
	  sum += temp2[a][h];
	}
      }
      */
      // ippsAdd
      for (int a = nangles - 2; a >= 0; a--) {
	for (int h = 0; h < nh; h++)
	  pixels[a + 1][h][v] = pixels[a][h][v] + temp1[a][h];
      }
    }
  }
}

void CCPi::remove_column_ring_artefacts(pixel_data &pixels,
					const sl_int nangles, const sl_int nh,
					const sl_int nv)
{
  //sl_int nh = get_num_h_pixels(); // nxi
  //sl_int nv = get_num_v_pixels();
  //sl_int nangles = get_num_angles(); // ny
  //pixel_data &pixels = get_pixel_data();
  const int crop_left = 0;
  const int crop_right = 0;
  const int nx = nh - (crop_left + crop_right);
  real_1d vec_res(nx);
  real_1d vec64(nx);
  pixel_1d vec32(nx);
  // ippsZero
  for (int v = 0; v < nv; v++) {
    for (int h = 0; h < nx; h++)
      vec_res[h] = real(0.0);
    // Todo - don't need vec64 really
    for (int a = 0; a < nangles; a++) {
      // ippsConvert
      for (int h = 0; h < nx; h++)
	vec64[h] = real(pixels[a][h + crop_left][v]);
      // ippsAdd
      for (int h = 0; h < nx; h++)
	vec_res[h] += vec64[h];
    }
    // ippsDivC
    for (int h = 0; h < nx; h++)
      vec_res[h] /= real(nangles);
    // ippsConvert
    for (int h = 0; h < nx; h++)
      vec32[h] = pixel_type(vec_res[h]);
    // ippsSub
    for (int a = 0; a < nangles; a++) {
      for (int h = 0; h < nx; h++)
	pixels[a][h + crop_left][v] -= vec32[h];
    }
  }
}

void CCPi::remove_aml_ring_artefacts(pixel_data &pixels, const sl_int nangles,
				     const sl_int nh, const sl_int nv,
				     const real param_n, const real param_r,
				     const int num_series)
{
  //sl_int nh = get_num_h_pixels(); // nxi
  //sl_int nv = get_num_v_pixels();
  //sl_int nangles = get_num_angles(); // ny
  if (param_n < -1e-10) {
    report_error("Ring artefacts - wrong param N");
    return;
  }
  if (param_r < 1e-8) {
    report_error("Ring artefacts - wrong param R");
    return;
  }
  if (num_series < 1 or num_series > 100) {
    report_error("High peaks - num series out of range");
    return;
  }
  //pixel_data &pixels = get_pixel_data();
  const int crop_left = 0;
  const int crop_right = 0;
  const int nx = nh - (crop_left + crop_right);
  // based on Applied Mathematics Letters v23 p1489
  // Init
  real_3d mata(boost::extents[num_series][nx][nx]);
  real_2d vecRS(boost::extents[2 * num_series][nangles]);
  {
    // Todo - can combine a lot of this
    for(int s = 1; s < num_series; s++) {
      // ippsVectorSlope
      real slope = real(s) * M_PI / real(nangles);
      for (int i = 0; i < nangles; i++)
	vecRS[0][i] = real(i) * slope; // + 0.0
      // ippsSinCos
      for (int i = 0; i < nangles; i++) {
	vecRS[2 * s - 1][i] = std::sin(vecRS[0][i]);
	vecRS[2 * s - 0][i] = std::cos(vecRS[0][i]);
      }
      // ippsMulC
      real sqrt2 = std::sqrt(real(2.0));
      for (int i = 0; i < nangles; i++)
	vecRS[2 * s - 1][i] *= sqrt2;
      for (int i = 0; i < nangles; i++)
	vecRS[2 * s - 0][i] *= sqrt2;
    }
    // ippsSet
    for (int i = 0; i < nangles; i++)
      vecRS[0][i] = 1.0;
    real_1d vec_exp(2 * nx + 2);
    for (int s = 1; s <= num_series; s++) {
      real alpha = (param_n + param_r) * real(s * s);
      real mdiv = std::sqrt(alpha * (alpha + real(4.0)));
      real vec1 = std::sqrt(alpha) / real(2.0);
      //ippsAsinh
      real vec2 = asinh(vec1); // ln(x + sqrt(x * x + 1))
      real tau = 2.0 * vec2;
      // ippsVectorSlope
      for (int i = 0; i < 2 * nx + 2; i++)
	vec_exp[i] = -tau * real(i); // + 0.0
      // ippsExp
      for (int i = 0; i < 2 * nx + 2; i++)
	vec_exp[i] = std::exp(vec_exp[i]);
      int ii;
      int jj;
      for (int i = 0; i < nx; i++) {
	for (int j = 0; j < nx; j++) {
	  if (j > i) {
	    ii = i;
	    jj = j;
	  } else {
	    ii = j;
	    jj = i;
	  }
	  ii++;
	  jj++;
	  mata[s - 1][i][j] = vec_exp[abs(jj - ii)] *
	    (1.0 + vec_exp[2 * nx - 2 * jj + 1]) * (1.0 + vec_exp[2 * ii - 1]);
	}
      }
      // ippsDivC
      real dv = mdiv * (1.0 - vec_exp[2 * nx]);
      for (int i = 0; i < nx; i++)
	for (int j = 0; j < nx; j++)
	  mata[s - 1][i][j] /= dv;
    }
  }
  // per slice code
  real_1d vec_res(nx);
  real_1d vec64(nx);
  pixel_1d vec32(nx);
  int n = 2 * num_series - 1;
  for (int v = 0; v < nv; v++) {
    for (int s = 0; s < n; s++) {
      int ss = (s + 2) / 2;
      for (int h = 0; h < nx; h++)
	vec_res[h] = real(0.0);
      for (int a = 0; a < nangles; a++) {
	// ippsConvert
	for (int h = 0; h < nx; h++)
	  vec64[h] = real(pixels[a][h + crop_left][v]);
	// ippsMul
	for (int h = 0; h < nx; h++)
	  vec64[h] *= vecRS[s][a];
	// ippsAdd
	for (int h = 0; h < nx; h++)
	  vec_res[h] += vec64[h];
      }
      // ippsDivC
      for (int h = 0; h < nx; h++)
	vec_res[h] /= real(nangles);
      // ippsMulC
      for (int h = 0; h < nx; h++)
	vec64[h] = param_n * real(ss * ss) * vec_res[h];
      // ippsAdd
      for (int h = 0; h < nx - 1; h++)
	vec64[h] += vec_res[h];
      // ippsSub
      for (int h = 0; h < nx - 1; h++)
	vec64[h] -= vec_res[h + 1];
      // ippsAdd
      for (int h = 1; h < nx; h++)
	vec64[h] += vec_res[h];
      // ippsSub
      for (int h = 0; h < nx - 1; h++)
	vec64[h + 1] -= vec_res[h];
      // ippsDotProd
      // Todo - we don't seem to use mata[>0][][], why?, and can we avoid calc?
      for (int i = 0; i < nx; i++) {
	real dp = 0.0;
	for (int h = 0; h < nx; h++)
	  dp += vec64[h] * mata[0][i][h];
	vec_res[i] = dp;
      }
      // ippsConvert
      for (int h = 0; h < nx; h++)
	vec32[h] = pixel_type(vec_res[h]);
      for (int a = 0; a < nangles; a++) {
	// ippsAddProductC
	for (int h = 0; h < nx; h++)
	  pixels[a][h + crop_left][v] -= pixel_type(vecRS[s][a]) * vec32[h];
      }
    }
  }
}

// Todo - virtual fns and Diamond sub-classes for each alg?
void CCPi::Diamond::ring_artefact_removal(const ring_artefact_alg alg,
					  const real param_n,
					  const real param_r,
					  const int num_series)
{
  // ported from ring_artefacts routine in Manchester/Diamond FBP code
  if (alg == ring_artefacts_column)
    remove_column_ring_artefacts(get_pixel_data(), get_num_angles(),
				 get_num_h_pixels(), get_num_v_pixels());
  else if (alg == ring_artefacts_aml)
    remove_aml_ring_artefacts(get_pixel_data(), get_num_angles(),
			      get_num_h_pixels(), get_num_v_pixels(),
			      param_n, param_r, num_series);
}
