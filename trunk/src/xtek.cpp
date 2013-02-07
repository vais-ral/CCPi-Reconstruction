
#include <iostream>
#include <omp.h>
#include "base_types.hpp"
#include "instruments.hpp"
#include "project_line.hpp"
#include "xtek_b.hpp"
#include "xtek_f.hpp"

// Nikon XTek instrument
// 360 degree clockwise sample rotations about vertical axis, cone beam

bool CCPi::Nikon_XTek::setup_experimental_geometry(const std::string file,
						   const bool phantom)
{
  if (phantom)
    return create_phantom();
  else
    return read_config_file(file);
}

bool CCPi::Nikon_XTek::create_phantom()
{
  source_x = -250.0;
  source_y = 0.0;
  source_z = 0.0;
  detector_x = 737.0;
  horizontal_pixels = new real[2000];
  int i;
  for (i = 0; i < 2000; i++) {
    real p = -116.5225 + i * 0.254;
    if (p <= 116.5225 + 0.001)
      horizontal_pixels[i] = p;
    else
      break;
  }
  n_horizontal_pixels = i;
  vertical_pixels = new real[2000];
  for (i = 0; i < 2000; i++) {
    real p = -92.3925 + i * 0.254;
    if (p <= 92.3925 + 0.001)
      vertical_pixels[i] = p;
    else
      break;
  }
  n_vertical_pixels = i;
  // 501 values from 0 to 2pi
  //geom.angles = linspace(0,2*pi,501);
  // lose the 2pi which is a duplicate of 0
  //geom.angles = geom.angles(1:500);
  angles = new real[500];
  real step = 2.0 * M_PI / 250;
  for (i = 0; i < 250; i++)
    angles[i] = i * step;
  n_angles = 250;
  return true;
}

bool CCPi::Nikon_XTek::read_config_file(const std::string file)
{
  std::cerr << "Todo\n";
  return false;
}

bool CCPi::Nikon_XTek::read_scans(const bool phantom)
{
  if (phantom)
    return build_phantom();
  else
    return read_images();
}

bool CCPi::Nikon_XTek::build_phantom()
{
  // build small voxel set to project to produce initial pixels
  int nx = 500;
  int ny = 500;
  int nz = 500;

  real mask_radius = -source_x
    * std::sin(std::atan(horizontal_pixels[n_horizontal_pixels - 1] /
			 (detector_x - source_x)));
  real voxel_size[3];
  voxel_size[0] = (2 * mask_radius / nx);
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

  long n_rays = n_angles * n_vertical_pixels * n_horizontal_pixels;
  pixel_type *pixels = new pixel_type[n_rays];
  // perform projection step
  forward_project(source_x, source_y, source_z, detector_x, horizontal_pixels,
		  vertical_pixels, angles, pixels, x, n_angles,
		  n_horizontal_pixels, n_vertical_pixels, image_offset,
		  voxel_size, nx, ny, nz);
  delete [] x;
  // Todo - could do with adding noise.
  return true;
}

bool CCPi::Nikon_XTek::read_images()
{
  std::cerr << "Todo\n";
  return false;
}

void CCPi::Nikon_XTek::forward_project(voxel_data &voxels,
				       const real origin[3],
				       const real width[3])
{
  const voxel_data::size_type *s = voxels.shape();
  forward_project(source_x, source_y, source_z, detector_x, horizontal_pixels,
		  vertical_pixels, angles, pixels, voxels.data(), n_angles,
		  n_horizontal_pixels, n_vertical_pixels, origin, width,
		  (int)s[0], (int)s[1], (int)s[2]);
}

void CCPi::Nikon_XTek::backward_project(voxel_data &voxels,
					const real origin[3],
					const real width[3])
{
  const voxel_data::size_type *s = voxels.shape();
  backward_project(source_x, source_y, source_z, detector_x, horizontal_pixels,
		   vertical_pixels, angles, pixels, voxels.data(), n_angles,
		   n_horizontal_pixels, n_vertical_pixels, origin, width,
		   (int)s[0], (int)s[1], (int)s[2]);
}
