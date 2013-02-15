
#include <fstream>
#include <cstring>
#include <cctype>
#include <omp.h>
#include "base_types.hpp"
#include "utils.hpp"
#include "instruments.hpp"
#include "project_line.hpp"
#include "xtek_b.hpp"
#include "xtek_f.hpp"
#include "tiffio.h"

// Nikon XTek instrument
// 360 degree clockwise sample rotations about vertical axis, cone beam

bool CCPi::Nikon_XTek::setup_experimental_geometry(const std::string path,
						   const std::string file,
						   const bool phantom)
{
  if (phantom)
    return create_phantom();
  else
    return read_config_file(path, file);
}

bool CCPi::Nikon_XTek::create_phantom()
{
  source_x = -250.0;
  source_y = 0.0;
  source_z = 0.0;
  detector_x = 737.0;
  horizontal_pixels = new real[1000];
  int i;
  for (i = 0; i < 1000; i++) {
    //real p = -100.0 + i * 0.250;
    real p = -100.0 + i * 0.400;
    if (p <= 100.0 + 0.001)
      horizontal_pixels[i] = p;
    else
      break;
  }
  n_horizontal_pixels = i;
  vertical_pixels = new real[1000];
  for (i = 0; i < 1000; i++) {
    //real p = -100.0 + i * 0.250;
    real p = -100.0 + i * 0.400;
    if (p <= 100.0 + 0.001)
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

bool CCPi::Nikon_XTek::read_config_file(const std::string path,
					const std::string file)
{
  bool ok = false;
  // read Xtek .xtekct file
  std::string data_file;
  combine_path_and_name(path, file, data_file);
  std::ifstream input(data_file.c_str());
  if (input.good()) {
    char line[256];
    real xpixel_size = 0.0;
    real ypixel_size = 0.0;
    real init_angle = 0.0;
    white_level = 0.0;
    input.getline(line, 256);
    if (strncmp(line, "[XTekCT]", 8) == 0) {
      ok = true;
      int ndata = 0;
      while (!input.eof()) {
	input.getline(line, 256);
	if (strncmp(line, "SrcToObject", 11) == 0) {
	  source_x = - std::atof(&line[12]);
	  source_y = 0.0;
	  source_z = 0.0;
	  ndata++;
	} else if (strncmp(line, "SrcToDetector", 13) == 0) {
	  detector_x = std::atof(&line[14]);
	  ndata++;
	} else if (strncmp(line, "MaskRadius", 10) == 0) {
	  mask_radius = std::atof(&line[11]);
	  ndata++;
	} else if (strncmp(line, "DetectorPixelsX", 15) == 0) {
	  n_horizontal_pixels = std::atoi(&line[16]);
	  ndata++;
	} else if (strncmp(line, "DetectorPixelsY", 15) == 0) {
	  n_vertical_pixels = std::atoi(&line[16]);
	  ndata++;
	} else if (strncmp(line, "DetectorPixelSizeX", 18) == 0) {
	  xpixel_size = std::atof(&line[19]);
	  ndata++;
	} else if (strncmp(line, "DetectorPixelSizeY", 18) == 0) {
	  ypixel_size = std::atof(&line[19]);
	  ndata++;
	} else if (strncmp(line, "Projections", 11) == 0) {
	  n_angles = std::atoi(&line[12]);
	  ndata++;
	} else if (strncmp(line, "InitialAngle", 12) == 0) {
	  init_angle = std::atof(&line[13]);
	  ndata++;
	} else if (strncmp(line, "WhiteLevel", 10) == 0) {
	  white_level = std::atof(&line[11]);
	  ndata++;
	} else if (strncmp(line, "Name", 4) == 0) {
	  while (isspace(line[strlen(line) - 1]))
	    line[strlen(line) - 1] = '\0';
	  basename = &line[5];
	  ndata++;
	}
	if (ndata == 11)
	  break;
      }
      if (ndata == 11) {
	detector_x += source_x;
	if (basename == "") {
	  ok = false;
	  std::cerr << "No name found for tiff files\n";
	} else if (n_horizontal_pixels > 0 and n_vertical_pixels > 0 and
		   n_angles > 0) {
	  horizontal_pixels = new real[n_horizontal_pixels];
	  vertical_pixels = new real[n_vertical_pixels];
	  angles = new real[n_angles];
	  real pixel_base = -((n_horizontal_pixels - 1) * xpixel_size / 2.0);
	  for (int i = 0; i < n_horizontal_pixels; i++)
	    horizontal_pixels[i] = pixel_base + i * xpixel_size;
	  pixel_base = -((n_vertical_pixels - 1) * ypixel_size / 2.0);
	  for (int i = 0; i < n_vertical_pixels; i++)
	    vertical_pixels[i] = pixel_base + i * ypixel_size;
	  std::string ctfile;
	  combine_path_and_name(path, "_ctdata.txt", ctfile);
	  ok = read_angles(ctfile, init_angle);
	} else {
	  ok = false;
	  std::cerr << "Negative pixel values\n";
	}
      } else {
	ok = false;
	std::cerr << "Failed to locate all geometry data in xtekct file\n";
      }
    } else
      std::cerr << "Incorrect header on xtekct file\n";
    input.close();
  } else
    std::cerr << "Open " << file << " failed\n";
  return ok;
}

bool CCPi::Nikon_XTek::read_angles(const std::string datafile,
				   const real init_angle)
{
  // since we read actual angles the answer is probably no
  //std::cout << "Should we use initial angle?\n";
  std::ifstream data(datafile.c_str());
  if (data.good()) {
    // sets of 3 items terminated by ^M, second is the angle.
    real tmp;
    for (int i = 0; i < n_angles; i++) {
      data >> tmp;
      data >> angles[i];
      data >> tmp;
    }
    data.close();
    return true;
  } else {
    std::cerr << "Error opening ctdata file\n";
    return false;
  }
}

bool CCPi::Nikon_XTek::finish_voxel_geometry(real voxel_origin[3],
					     real voxel_size[3],
					     const voxel_data &voxels) const
{
  const voxel_data::size_type *s = voxels.shape();
  real size = 2.0 * mask_radius / real(s[0]);
  voxel_size[0] = size;
  voxel_size[1] = size;
  voxel_size[2] = size;
  voxel_origin[0] = -voxel_size[0] * real(s[0]) / 2.0 + offset[0];
  voxel_origin[1] = -voxel_size[1] * real(s[1]) / 2.0 + offset[1];
  voxel_origin[2] = -voxel_size[2] * real(s[2]) / 2.0 + offset[2];
  return true;
}

bool CCPi::Nikon_XTek::read_scans(const std::string path, const bool phantom)
{
  if (phantom)
    return build_phantom();
  else
    return read_images(path);
}

bool CCPi::Nikon_XTek::build_phantom()
{
  // build small voxel set to project to produce initial pixels
  int nx = 500;
  int ny = 500;
  int nz = 500;

  mask_radius = -source_x
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
  pixel_data = new pixel_type[n_rays];
  // perform projection step
  forward_project(source_x, source_y, source_z, detector_x, horizontal_pixels,
		  vertical_pixels, angles, pixel_data, x, n_angles,
		  n_horizontal_pixels, n_vertical_pixels, image_offset,
		  voxel_size, nx, ny, nz);
  delete [] x;
  // Todo - could do with adding noise.
  offset[0] = 0.0;
  offset[1] = 0.0;
  offset[2] = 0.0;
  return true;
}

bool CCPi::Nikon_XTek::read_images(const std::string path)
{
  bool ok = true;
  long n_rays = n_angles * n_vertical_pixels * n_horizontal_pixels;
  pixel_data = new pixel_type[n_rays];
  std::string pathbase;
  combine_path_and_name(path, basename, pathbase);
  char index[8];
  for (int i = 0; i < n_angles; i++) {
    std::cout << "loop " << i << '\n';
    snprintf(index, 8, "%04d", i + 1);
    std::string name = pathbase + index + ".tif";
    TIFF *tif = TIFFOpen(name.c_str(), "r");
    if (tif == 0) {
      ok = false;
      std::cerr << "Error opening " << name << '\n';
    } else {
      if (TIFFIsTiled(tif)) {
	ok = false;
        std::cerr << "tiled image not supported in XTek reader\n";
      } else {
	uint16 bps;
	if (TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bps) == 0) {
	  std::cerr << "TIFF error reading bits per sample\n";
	  ok = false;
	} else if (bps != 16) {
	  std::cerr << "TIFF is not 16bit data\n";
	  ok = false;
	} else {
	  uint16 photo;
	  if (TIFFGetField(tif, TIFFTAG_PHOTOMETRIC, &photo) == 0) {
	    std::cerr << "Error getting TIFF photometric info\n";
	    ok = false;
	  } else if (photo != PHOTOMETRIC_MINISBLACK) {
	    std::cerr << "TIFF photometric type not supported by XTek reader\n";
	    ok = false;
	  } else {
	    if (TIFFIsMSB2LSB(tif) == 0) {
	      std::cerr << "TIFF is not MSB to LSB\n";
	      ok = false;
	    } else {
	      // I'm assuming orientation would be 1, (0,0) in top left corner
	      // but XTek tiffs don't usually work with TIFFTAG_ORIENTATION
	      uint32 w;
	      if (TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w) == 0) {
		std::cerr << "Error getting TIFF width\n";
		ok = false;
	      } else if ((int)w != n_horizontal_pixels) {
		std::cerr << "Image width mismatch\n";
		ok = false;
	      } else {
		tsize_t strip_size = TIFFStripSize(tif);
		int nstrips = TIFFNumberOfStrips(tif);
		char *buf = (char *)_TIFFmalloc(strip_size * nstrips);
		int offset = 0;
		int result;
		for (int count = 0; (count < nstrips and ok); count++){
		  if ((result = TIFFReadEncodedStrip(tif, count, buf + offset,
						     strip_size)) == -1) {
		    std::cerr << "Read error in " << name << '\n';
		    ok = false;
		  } else
		    offset += result;
		}
		int extra = offset % 2;
		offset /= 2;
		if (extra == 1) {
		  std::cerr << "Odd number of bytes for 16 bit image\n";
		  ok = false;
		} else if (offset != n_horizontal_pixels * n_vertical_pixels) {
		  std::cerr << "Image size mismatch\n";
		  ok = false;
		} else {
		  // (0,0) at top left means vertical order needs to be reversed
		  // for voxel box where the origin is at the bottom
		  // It looks like strips are [h][v] col-major
		  // storage order is [horiz][vert][angles] fortran order
		  // for compatibility with Matlab.
		  uint16 *b = (uint16 *)buf;
		  long angle_offset = i * n_horizontal_pixels
		    * n_vertical_pixels;
		  for (int v = 0; v < n_vertical_pixels; v++) {
		    for (int h = 0; h < n_horizontal_pixels; h++) {
		      pixel_data[angle_offset] =
			real(b[(n_vertical_pixels - v - 1)
			       * n_horizontal_pixels + h]);
		      angle_offset++;
		    }
		  }
		}
		_TIFFfree(buf);
	      }
	    }
	  }
	}
      }
      TIFFClose(tif);
    }
  }
  if (ok) {
    /*
    real max_v = 0.0;
    for (long j = 0; j < n_rays; j++)
      if (max_v < pixel_data[j])
	max_v = pixel_data[j];
    if (max_v > white_level + 0.01)
      std::cout << "Values exceed white level\n";
    else
      max_v = white_level;
    */
    real max_v = 65535.0;
    // scale and take -ve log, due to exponential extinction in sample.
    for (long j = 0; j < n_rays; j++)
      pixel_data[j] = - std::log(pixel_data[j] / max_v);
  }
  return ok;
}

void CCPi::Nikon_XTek::apply_beam_hardening()
{
  long n_rays = n_angles * n_vertical_pixels * n_horizontal_pixels;
  for (long i = 0; i < n_rays; i++)
    pixel_data[i] = pixel_data[i] * pixel_data[i];
}

long CCPi::Nikon_XTek::get_data_size() const
{
  return n_angles * n_vertical_pixels * n_horizontal_pixels;
}

pixel_type *const CCPi::Nikon_XTek::get_pixel_data() const
{
  return pixel_data;
}

void CCPi::Nikon_XTek::forward_project(pixel_type *pixels,
				       voxel_type *const voxels,
				       const real origin[3],
				       const real width[3], const int nx,
				       const int ny, const int nz) const
{
  forward_project(source_x, source_y, source_z, detector_x, horizontal_pixels,
		  vertical_pixels, angles, pixels, voxels, n_angles,
		  n_horizontal_pixels, n_vertical_pixels, origin, width,
		  nx, ny, nz);
}

void CCPi::Nikon_XTek::backward_project(pixel_type *pixels,
					voxel_type *const voxels,
					const real origin[3],
					const real width[3], const int nx,
					const int ny, const int nz) const
{
  backward_project(source_x, source_y, source_z, detector_x, horizontal_pixels,
		   vertical_pixels, angles, pixels, voxels, n_angles,
		   n_horizontal_pixels, n_vertical_pixels, origin, width,
		   nx, ny, nz);
}

void CCPi::Nikon_XTek::backward_project(voxel_type *const voxels,
					const real origin[3],
					const real width[3], const int nx,
					const int ny, const int nz) const
{
  backward_project(source_x, source_y, source_z, detector_x, horizontal_pixels,
		   vertical_pixels, angles, pixel_data, voxels, n_angles,
		   n_horizontal_pixels, n_vertical_pixels, origin, width,
		   nx, ny, nz);
}
