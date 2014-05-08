
#include <iostream>
#include "base_types.hpp"
#include "instruments.hpp"
#include "nexus.hpp"
#include "utils.hpp"

bool CCPi::Diamond::setup_experimental_geometry(const std::string path,
						const std::string file,
						const bool phantom)
{
  if (phantom)
    return create_phantom();
  else {
    name = file;
    return true;
  }
}

bool CCPi::Diamond::create_phantom()
{
  //set_source(-250.0, 0.0, 0.0);
  //set_detector(737.0);
  int c;
  for (c = 0; c < 1000; c++) {
    //real p = -100.0 + c * 0.250;
    real p = -99.8 + real(c) * 0.400;
    if (p >= 100.0 + 0.001)
      break;
  }
  real_1d &h_pixels = set_h_pixels(c);
  for (int i = 0; i < c; i++) {
    //real p = -100.0 + i * 0.250;
    h_pixels[i] = -99.8 + real(i) * 0.400;
  }
  for (c = 0; c < 1000; c++) {
    //real p = -100.0 + real(c) * 0.250;
    real p = -99.8 + real(c) * 0.400;
    if (p >= 100.0 + 0.001)
      break;
  }
  real_1d &v_pixels = set_v_pixels(c);
  for (int i = 0; i < c; i++) {
    //real p = -100.0 + real(i) * 0.250;
    v_pixels[i] = -99.8 + real(i) * 0.400;
  }
  // 501 values from 0 to 2pi
  //geom.angles = linspace(0,2*pi,501);
  // lose the 2pi which is a duplicate of 0
  //geom.angles = geom.angles(1:500);
  const int nangles = 250;
  real_1d &pangles = set_phi(nangles);
  real step = 2.0 * M_PI / real(nangles);
  for (int i = 0; i < nangles; i++)
    pangles[i] = real(i) * step;
  return true;
}

bool CCPi::Diamond::read_data_size(const std::string path,
				   const bool phantom)
{
  // phantom already done by setup
  if (phantom)
    return true;
  else {
    bool ok = false;
    std::string fullname;
    combine_path_and_name(path, name,fullname);
    pixel_type *pixels = 0;
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
      // store data in class
      real_1d &h_pixels = set_h_pixels(nh_pixels);
      h_pixels[0] = - ((nh_pixels - 1) * hsize) / real(2.0);
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
  voxel_data x(boost::extents[nx][ny][nz], boost::fortran_storage_order());
  //sl_int n_vox = nx * ny * nz;
  for (sl_int i = 0; i < nz; i++)
    for (sl_int j = 0; j < ny; j++)
      for (sl_int k = 0; k < nx; k++)
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
  pixel_type *pixels = create_pixel_data();
  // perform projection step
  forward_project(pixels, x, image_offset, voxel_size, nx, ny, nz);
  //delete [] x;
  return true;
}

bool CCPi::Diamond::read_data(const std::string path, const int offset,
			      const int block_size, const bool first)
{
  bool ok = false;
  std::string fullname;
  combine_path_and_name(path, name,fullname);
  pixel_type *pixels = 0;
  int nh_pixels = 0;
  int nv_pixels = 0;
  std::vector<real> angles(2 * get_num_angles());
  int nangles = 0;
  real hsize = 0.0;
  real vsize = 0.0;
  int nv = get_num_v_pixels();
  int nh = get_num_h_pixels();
  sl_int sz =  sl_int(nv) * sl_int(nh);
  if (first)
    pixels = new pixel_type[get_num_angles() * sz];
  else
    pixels = get_pixel_data();
  pixel_2d i_dark(boost::extents[nv][nh]);
  pixel_2d f_dark(boost::extents[nv][nh]);
  pixel_2d i_bright(boost::extents[nv][nh]);
  pixel_2d f_bright(boost::extents[nv][nh]);
  for (sl_int i = 0; i < nv; i++) {
    for (sl_int j = 0; j < nh; j++) {
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
    sl_int n_rays = sl_int(nangles) * sz;
    if (first)
      set_pixel_data(pixels, n_rays);
    const real_1d &angles = get_phi();
    // linear interpolate bright/dark frames. Todo - something else?
    // Todo, also use initial final bright/dark angles? rather than assuming
    // initial/final sample angles?
    sl_int n = sz;
    pixel_2d dark(boost::extents[nv][nh]);
    pixel_2d bright(boost::extents[nv][nh]);
    for (int i = 0; i < nangles; i++) {
      // Based on fbp code, interpolate bright/dark
      // w = (angles[i] - angles[0]) / (angles[nangles - 1] - angles[0])?
      real w = angles[i] / angles[nangles - 1];
      for (sl_int j = 0; j < nv; j++)
	for (sl_int k = 0; k < nh; k++)
	  dark[j][k] = i_dark[j][k] * (real(1.0) - w) + f_dark[j][k] * w;
      for (sl_int j = 0; j < nv; j++)
	for (sl_int k = 0; k < nh; k++)
	  bright[j][k] = i_bright[j][k] * (real(1.0) - w) + f_bright[j][k] * w;
      // subtract dark from data/bright
      // and clamp min data/bright value to 0.1
      for (sl_int j = 0; j < nv; j++) {
	for (sl_int k = 0; k < nh; k++) {
	  bright[j][k] -= dark[j][k];
	  if (bright[j][k] < real(0.1))
	    bright[j][k] = 0.1;
	}
      }
      for (sl_int j = 0; j < nv; j++) {
	for (sl_int k = 0; k < nh; k++) {
	  pixels[k + j * nh + i * n] -= dark[j][k];
	  if (pixels[k + j * nh + i * n] < real(0.1))
	    pixels[k + j * nh + i * n] = 0.1;
	}
      }
      // scale each data pixel by bright pixel
      for (sl_int j = 0; j < nv; j++)
	for (sl_int k = 0; k < nh; k++)
	  pixels[k + j * nh + i * n] /= bright[j][k];
    }
    // take -ve log, due to exponential extinction in sample.
    for (sl_int j = 0; j < n_rays; j++)
      pixels[j] = - std::log(pixels[j]);
    //find_centre(get_num_v_pixels() / 2 + 1);
  }
  return ok;
}

bool CCPi::Diamond::finish_voxel_geometry(real voxel_origin[3],
					  real voxel_size[3], const int nx,
					  const int ny, const int nz) const
{
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
  voxel_size[0] = hsize;
  voxel_size[1] = hsize;
  voxel_size[2] = vsize;
  voxel_origin[0] = -voxel_size[0] * real(nx) / real(2.0); // + offset[0];
  voxel_origin[1] = -voxel_size[1] * real(ny) / real(2.0); // + offset[1];
  voxel_origin[2] = -voxel_size[2] * real(nz) / real(2.0); // + offset[2];
  return true;
}

void CCPi::Diamond::apply_beam_hardening()
{
  // Todo - does this belong in the base class?
  sl_int n_rays = get_data_size();
  pixel_type *pixels = get_pixel_data();
  for (sl_int i = 0; i < n_rays; i++)
    pixels[i] = pixels[i] * pixels[i];
}
