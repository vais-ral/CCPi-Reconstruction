
#include <fstream>
#include <cstring>
#include <cctype>
#include <omp.h>
#include "base_types.hpp"
#include "utils.hpp"
#include "instruments.hpp"
#include "project_line.hpp"
#include "cone_b.hpp"
#include "cone_f.hpp"
#include "tiff.hpp"
#include "timer.hpp"

#ifndef USE_TIMER
#  define USE_TIMER false
#endif // USE_TIMER

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
  set_source(-250.0, 0.0, 0.0);
  set_detector(737.0);
  real *h_pixels = new real[1000];
  int i;
  for (i = 0; i < 1000; i++) {
    //real p = -100.0 + i * 0.250;
    real p = -100.0 + i * 0.400;
    if (p <= 100.0 + 0.001)
      h_pixels[i] = p;
    else
      break;
  }
  set_h_pixels(h_pixels, i);
  real *v_pixels = new real[1000];
  for (i = 0; i < 1000; i++) {
    //real p = -100.0 + i * 0.250;
    real p = -100.0 + i * 0.400;
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
    real det_x = 0.0;
    real xpixel_size = 0.0;
    real ypixel_size = 0.0;
    real init_angle = 0.0;
    white_level = 0.0;
    int n_h_pixels = 0;
    int n_v_pixels = 0;
    int n_phi = 0;
    input.getline(line, 256);
    if (strncmp(line, "[XTekCT]", 8) == 0) {
      ok = true;
      int ndata = 0;
      while (!input.eof()) {
	input.getline(line, 256);
	if (strncmp(line, "SrcToObject", 11) == 0) {
	  set_source(- std::atof(&line[12]), 0.0, 0.0);
	  ndata++;
	} else if (strncmp(line, "SrcToDetector", 13) == 0) {
	  det_x = std::atof(&line[14]);
	  ndata++;
	} else if (strncmp(line, "MaskRadius", 10) == 0) {
	  mask_radius = std::atof(&line[11]);
	  ndata++;
	} else if (strncmp(line, "DetectorPixelsX", 15) == 0) {
	  n_h_pixels = std::atoi(&line[16]);
	  ndata++;
	} else if (strncmp(line, "DetectorPixelsY", 15) == 0) {
	  n_v_pixels = std::atoi(&line[16]);
	  ndata++;
	} else if (strncmp(line, "DetectorPixelSizeX", 18) == 0) {
	  xpixel_size = std::atof(&line[19]);
	  ndata++;
	} else if (strncmp(line, "DetectorPixelSizeY", 18) == 0) {
	  ypixel_size = std::atof(&line[19]);
	  ndata++;
	} else if (strncmp(line, "Projections", 11) == 0) {
	  n_phi = std::atoi(&line[12]);
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
	set_detector(det_x + get_source_x());
	if (basename == "") {
	  ok = false;
	  std::cerr << "No name found for tiff files\n";
	} else if (n_h_pixels > 0 and n_v_pixels > 0 and n_phi > 0) {
	  real *h_pixels = new real[n_h_pixels];
	  real *v_pixels = new real[n_v_pixels];
	  real pixel_base = -((n_h_pixels - 1) * xpixel_size / 2.0);
	  for (int i = 0; i < n_h_pixels; i++)
	    h_pixels[i] = pixel_base + i * xpixel_size;
	  set_h_pixels(h_pixels, n_h_pixels);
	  pixel_base = -((n_v_pixels - 1) * ypixel_size / 2.0);
	  for (int i = 0; i < n_v_pixels; i++)
	    v_pixels[i] = pixel_base + i * ypixel_size;
	  set_v_pixels(v_pixels, n_v_pixels);
	  std::string ctfile;
	  combine_path_and_name(path, "_ctdata.txt", ctfile);
	  ok = read_angles(ctfile, init_angle, n_phi);
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
				   const real init_angle, const int n)
{
  // since we read actual angles the answer is probably no
  //std::cout << "Should we use initial angle?\n";
  std::ifstream data(datafile.c_str());
  if (data.good()) {
    real *p = new real[n];
    // sets of 3 items terminated by ^M, second is the angle.
    real tmp;
    for (int i = 0; i < n; i++) {
      data >> tmp;
      data >> p[i];
      data >> tmp;
    }
    data.close();
    set_phi(p, n);
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

  real *h_pixels = get_h_pixels();
  int n_h_pixels = get_num_h_pixels();
  mask_radius = -get_source_x()
    * std::sin(std::atan(h_pixels[n_h_pixels - 1] /
			 (get_detector_x() - get_source_x())));
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

  pixel_type *pixels = create_pixel_data();
  // perform projection step
  instrument::forward_project(get_source_x(), get_source_y(), get_source_z(),
			      get_detector_x(), get_h_pixels(), get_v_pixels(),
			      get_phi(), get_theta(), pixels, x,
			      get_num_angles(), n_h_pixels, get_num_v_pixels(),
			      image_offset, voxel_size, nx, ny, nz);
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
  pixel_type *pixels = create_pixel_data();
  long n_rays = get_data_size();
  std::string pathbase;
  combine_path_and_name(path, basename, pathbase);
  char index[8];
  for (int i = 0; (i < get_num_angles() and ok); i++) {
    snprintf(index, 8, "%04d", i + 1);
    std::string name = pathbase + index + ".tif";
    long angle_offset = i * get_num_h_pixels() * get_num_v_pixels();
    ok = read_tiff(name, &pixels[angle_offset], get_num_h_pixels(),
		   get_num_v_pixels());
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
      pixels[j] = - std::log(pixels[j] / max_v);
  }
  return ok;
}

void CCPi::Nikon_XTek::apply_beam_hardening()
{
  // Todo - does this belong in the base class?
  long n_rays = get_data_size();
  pixel_type *pixels = get_pixel_data();
  for (long i = 0; i < n_rays; i++)
    pixels[i] = pixels[i] * pixels[i];
}

void CCPi::cone_beam::forward_project(pixel_type *pixels,
				      voxel_type *const voxels,
				      const real origin[3],
				      const real width[3], const int nx,
				      const int ny, const int nz) const
{
  timer fptime(USE_TIMER);
  instrument::forward_project(source_x, source_y, source_z, detector_x,
			      get_h_pixels(), get_v_pixels(), get_phi(),
			      get_theta(), pixels, voxels, get_num_angles(),
			      get_num_h_pixels(), get_num_v_pixels(), origin,
			      width, nx, ny, nz);
  fptime.accumulate();
  fptime.output(" forward projection");
}

void CCPi::cone_beam::backward_project(pixel_type *pixels,
				       voxel_type *const voxels,
				       const real origin[3],
				       const real width[3], const int nx,
				       const int ny, const int nz) const
{
  timer bptime(USE_TIMER);
  instrument::backward_project(source_x, source_y, source_z, detector_x,
			       get_h_pixels(), get_v_pixels(), get_phi(),
			       get_theta(), pixels, voxels, get_num_angles(),
			       get_num_h_pixels(), get_num_v_pixels(), origin,
			       width, nx, ny, nz);
  bptime.accumulate();
  bptime.output("backward projection");
}

void CCPi::cone_beam::backward_project(voxel_type *const voxels,
				       const real origin[3],
				       const real width[3], const int nx,
				       const int ny, const int nz) const
{
  timer bptime(USE_TIMER);
  instrument::backward_project(source_x, source_y, source_z, detector_x,
			       get_h_pixels(), get_v_pixels(), get_phi(),
			       get_theta(), get_pixel_data(), voxels,
			       get_num_angles(), get_num_h_pixels(),
			       get_num_v_pixels(), origin, width, nx, ny, nz);
  bptime.accumulate();
  bptime.output("backward projection");
}
