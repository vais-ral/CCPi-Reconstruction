
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
  real *h_pixels = new real[1000];
  int i;
  for (i = 0; i < 1000; i++) {
    //real p = -100.0 + i * 0.250;
    real p = -99.8 + i * 0.400;
    if (p <= 100.0 + 0.001)
      h_pixels[i] = p;
    else
      break;
  }
  set_h_pixels(h_pixels, i);
  real *v_pixels = new real[1000];
  for (i = 0; i < 1000; i++) {
    //real p = -100.0 + i * 0.250;
    real p = -99.8 + i * 0.400;
    if (p <= 100.0 + 0.001)
      v_pixels[i] = p;
    else
      break;
  }
  set_v_pixels(v_pixels, i);
  // 501 values from 0 to 2pi
  //geom.angles = linspace(0,2*pi,501);
  // lose the 2pi which is a duplicate of 0
  //geom.angles = geom.angles(1:500);
  real *pangles = new real[500];
  real step = 2.0 * M_PI / 250;
  for (i = 0; i < 250; i++)
    pangles[i] = i * step;
  set_phi(pangles, 250);
  return true;
}

bool CCPi::Diamond::read_scans(const std::string path, const bool phantom)
{
  if (phantom)
    return build_phantom();
  else
    return read_data(path);
}

bool CCPi::Diamond::build_phantom()
{
  // build small voxel set to project to produce initial pixels
  int nx = 500;
  int ny = 500;
  int nz = 500;

  real *h_pixels = get_h_pixels();
  real *v_pixels = get_v_pixels();
  int nh = get_num_h_pixels();
  int nv = get_num_v_pixels();
  real hmin = std::abs(h_pixels[0]);
  real hmax = h_pixels[nh - 1] + (h_pixels[1] - h_pixels[0]);
  real vmin = std::abs(v_pixels[0]);
  real vmax = v_pixels[nv - 1] + (v_pixels[1] - v_pixels[0]);
  real hlim = std::max(hmax, hmin);
  real vlim = std::max(vmax, vmin);
  real lim = std::max(hlim, vlim);
  real voxel_size[3];
  voxel_size[0] = 2.0 * lim / real(nx);
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
  long n_vox = nx * ny * nz;
  voxel_type *x = new voxel_type[n_vox];
  for (long i = 0; i < n_vox; i++)
    x[i] = 0.0;

  // add cubes - column major
  for (int i = 108-1; i < 189; i++) {
    for (int j = 108-1; j < 189; j++) {
      for (int k = 58-1; k < 139; k++) {
        x[k * nx * ny + j * nx + i] = 1;
      }
    }
  }
  for (int i = 190-1; i < 271; i++) {
    for (int j = 190-1; j < 271; j++) {
      for (int k = 140-1; k < 221; k++) {
        x[k * nx * ny + j * nx + i] = 1;
      }
    }
  }
  for (int i = 272-1; i < 353; i++) {
    for (int j = 272-1; j < 353; j++) {
      for (int k = 222-1; k < 303; k++) {
        x[k * nx * ny + j * nx + i] = 1;
      }
    }
  }

  pixel_type *pixels = create_pixel_data();
  // perform projection step
  forward_project(pixels, x, image_offset, voxel_size, nx, ny, nz);
  delete [] x;
  return true;
}

bool CCPi::Diamond::read_data(const std::string path)
{
  bool ok = true;
  std::string fullname;
  combine_path_and_name(path, name,fullname);
  pixel_type *pixels = 0;
  int nh_pixels = 0;
  int nv_pixels = 0;
  real *angles = 0;
  int nangles = 0;
  real hsize = 0.0;
  real vsize = 0.0;
  pixel_type *i_dark = 0;
  pixel_type *f_dark = 0;
  pixel_type *i_bright = 0;
  pixel_type *f_bright = 0;
  ok = read_NeXus(pixels, i_dark, f_dark, i_bright, f_bright, nh_pixels,
		  nv_pixels, angles, nangles, hsize, vsize, fullname, false);
  // store data in class
  real *h_pixels = new real[nh_pixels];
  h_pixels[0] = - ((nh_pixels - 1) * hsize) / 2.0;
  for (int i = 1; i < nh_pixels; i++)
    h_pixels[i] = h_pixels[0] + real(i) * hsize;
  set_h_pixels(h_pixels, nh_pixels);
  real *v_pixels = new real[nv_pixels];
  v_pixels[0] = - ((nv_pixels - 1) * vsize) / 2.0;
  for (int i = 1; i < nv_pixels; i++)
    v_pixels[i] = v_pixels[0] + real(i) * vsize;
  set_v_pixels(v_pixels, nv_pixels);
  set_phi(angles, nangles);
  long n_rays = long(nangles) * long(nh_pixels) * long(nv_pixels);
  set_pixel_data(pixels, n_rays);
  if (ok) {
    // linear interpolate bright/dark frames. Todo - something else?
    // Todo, also use initial final bright/dark angles? rather than assuming
    // initial/final sample angles?
    long n = nh_pixels * nv_pixels;
    pixel_type *dark = new pixel_type[n];
    pixel_type *bright = new pixel_type[n];
    for (int i = 0; i < nangles; i++) {
      // Based on fbp code, interpolate bright/dark
      // w = (angles[i] - angles[0]) / (angles[nangles - 1] - angles[0])?
      real w = angles[i] / angles[nangles - 1];
      for (long j = 0; j < n; j++)
	dark[j] = i_dark[j] * (1.0 - w) + f_dark[j] * w;
      for (long j = 0; j < n; j++)
	bright[j] = i_bright[j] * (1.0 - w) + f_bright[j] * w;
      // subtract dark from data/bright
      // and clamp min data/bright value to 0.1
      for (long j = 0; j < n; j++) {
	bright[j] -= dark[j];
	if (bright[j] < 0.1)
	  bright[j] = 0.1;
      }
      for (long j = 0; j < n; j++) {
	pixels[j + i * n] -= dark[j];
	if (pixels[j + i * n] < 0.1)
	  pixels[j + i * n] = 0.1;
      }
      // scale each data pixel by bright pixel
      for (long j = 0; j < n; j++)
	pixels[j + i * n] /= bright[j];
    }
    delete [] bright;
    delete [] dark;
    // take -ve log, due to exponential extinction in sample.
    for (long j = 0; j < n_rays; j++)
      pixels[j] = - std::log(pixels[j]);
    //find_centre(get_num_v_pixels() / 2 + 1);
  }
  delete [] f_bright;
  delete [] i_bright;
  delete [] f_dark;
  delete [] i_dark;
  return ok;
}

bool CCPi::Diamond::finish_voxel_geometry(real voxel_origin[3],
					  real voxel_size[3],
					  const voxel_data &voxels) const
{
  const voxel_data::size_type *s = voxels.shape();
  real *h_pixels = get_h_pixels();
  real *v_pixels = get_v_pixels();
  int nh = get_num_h_pixels();
  int nv = get_num_v_pixels();
  real hrange = std::abs(h_pixels[nh - 1] - h_pixels[0]);
  real hsize = hrange / (nh - 1);
  real vrange = std::abs(v_pixels[nv - 1] - v_pixels[0]);
  real vsize = vrange / (nv - 1);
  int hs = std::min(int(s[0]), int(s[1]));
  hsize *= (real(nh) / real(hs));
  vsize *= (real(nv) / real(s[2]));
  voxel_size[0] = hsize;
  voxel_size[1] = hsize;
  voxel_size[2] = vsize;
  voxel_origin[0] = -voxel_size[0] * real(s[0]) / 2.0; // + offset[0];
  voxel_origin[1] = -voxel_size[1] * real(s[1]) / 2.0; // + offset[1];
  voxel_origin[2] = -voxel_size[2] * real(s[2]) / 2.0; // + offset[2];
  return true;
}

void CCPi::Diamond::apply_beam_hardening()
{
  // Todo - does this belong in the base class?
  long n_rays = get_data_size();
  pixel_type *pixels = get_pixel_data();
  for (long i = 0; i < n_rays; i++)
    pixels[i] = pixels[i] * pixels[i];
}