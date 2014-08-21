
#include <iostream>
#include "base_types.hpp"
#include "fbp.hpp"
#include "instruments.hpp"
#include "nexus.hpp"
#include "utils.hpp"

bool CCPi::Diamond::setup_experimental_geometry(const std::string path,
						const std::string file,
						const real rotation_centre,
						const bool phantom)
{
  if (phantom)
    return create_phantom();
  else {
    name = file;
    return read_data_size(path, rotation_centre);
  }
}

bool CCPi::Diamond::create_phantom()
{
  //set_source(-250.0, 0.0, 0.0);
  //set_detector(737.0);
  int c;
  for (c = 0; c < 1000; c++) {
    //real p = -100.0 + c * 0.250;
    real p = -99.6 + real(c) * 0.800;
    if (p >= 100.0 + 0.001)
      break;
  }
  real_1d &h_pixels = set_h_pixels(c);
  for (int i = 0; i < c; i++) {
    //real p = -100.0 + i * 0.250;
    h_pixels[i] = -99.6 + real(i) * 0.800;
  }
  for (c = 0; c < 1000; c++) {
    //real p = -100.0 + real(c) * 0.250;
    real p = -99.6 + real(c) * 0.800;
    if (p >= 100.0 + 0.001)
      break;
  }
  real_1d &v_pixels = set_v_pixels(c);
  for (int i = 0; i < c; i++) {
    //real p = -100.0 + real(i) * 0.250;
    v_pixels[i] = -99.6 + real(i) * 0.800;
  }
  // 501 values from 0 to 2pi
  //geom.angles = linspace(0,2*pi,501);
  // lose the 2pi which is a duplicate of 0
  //geom.angles = geom.angles(1:500);
  const int nangles = 1000;
  real_1d &pangles = set_phi(nangles);
  real step = 2.0 * M_PI / real(nangles);
  for (int i = 0; i < nangles; i++)
    pangles[i] = real(i) * step;
  return true;
}

bool CCPi::Diamond::read_data_size(const std::string path,
				   const real rotation_centre)
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
		  false, false, 0, 10000);
#endif // HAS_NEXUS
  if (ok) {
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
  int nx = 500;
  int ny = 500;
  int nz = 500;

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
  int nh_pixels = 0;
  int nv_pixels = 0;
  std::vector<real> angles(2 * get_num_angles());
  int nangles = 0;
  real hsize = 0.0;
  real vsize = 0.0;
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
#ifdef HAS_NEXUS
  ok = read_NeXus(pixels, i_dark, f_dark, i_bright, f_bright, nh_pixels,
		  nv_pixels, angles, nangles, hsize, vsize, fullname,
		  false, true, offset, block_size);
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
      // w = (angles[i] - angles[0]) / (angles[nangles - 1] - angles[0])?
      real w = angles[i] / angles[nangles - 1];
      for (sl_int k = 0; k < nh; k++)
	for (sl_int j = 0; j < nv; j++)
	  dark[k][j] = i_dark[k][j] * (real(1.0) - w) + f_dark[k][j] * w;
      for (sl_int k = 0; k < nh; k++)
	for (sl_int j = 0; j < nv; j++)
	  bright[k][j] = i_bright[k][j] * (real(1.0) - w) + f_bright[k][j] * w;
      // subtract dark from data/bright
      // and clamp min data/bright value to 0.1
      for (sl_int k = 0; k < nh; k++) {
	for (sl_int j = 0; j < nv; j++) {
	  bright[k][j] -= dark[k][j];
	  if (bright[k][j] < real(0.1))
	    bright[k][j] = 0.1;
	}
      }
      for (sl_int k = 0; k < nh; k++) {
	for (sl_int j = 0; j < nv; j++) {
	  pixels[i][k][j] -= dark[k][j];
	  if (pixels[i][k][j] < real(0.1))
	    pixels[i][k][j] = 0.1;
	}
      }
      // scale each data pixel by bright pixel
      for (sl_int k = 0; k < nh; k++)
	for (sl_int j = 0; j < nv; j++)
	  pixels[i][k][j] /= bright[k][j];
    }
    // take -ve log, due to exponential extinction in sample.
    for (int i = 0; i < nangles; i++)
      for (sl_int k = 0; k < nh; k++)
	for (sl_int j = 0; j < nv; j++)
	  pixels[i][k][j] = - std::log(pixels[i][k][j]);
    //find_centre(get_num_v_pixels() / 2 + 1);
  }
  return ok;
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
