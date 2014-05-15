
#include <iostream>
#include <fstream>
#include <cstring>
#include <cctype>
#include "base_types.hpp"
#include "utils.hpp"
#include "instruments.hpp"
#include "tiff.hpp"
#include "ui_calls.hpp"

// Nikon XTek instrument
// 360 degree clockwise sample rotations about vertical axis, cone beam

static pixel_type linear(const real x1, const real x2, const pixel_type f1,
			 const pixel_type f2, const real x);
static pixel_type bilinear(const real x1, const real x2, const real y1,
			   const real y2, const pixel_type f11,
			   const pixel_type f12, const pixel_type f21,
			   const pixel_type f22, const real x, const real y);
static pixel_type angles_linear(const int ph1, const real_1d &h,
				const int v_slice, const real_1d &angles,
				const pixel_data &data, const real new_h,
				const real new_angle, const int nh,
				const int na, const int nv);
static pixel_type angles_bilinear(const int ph1, const real_1d &h,
				  const int v_slice, const real_1d &angles,
				  const pixel_data &data, const real new_h,
				  const real new_angle, const int nh,
				  const int na, const int nv);
static pixel_type interpolate2D(const real_1d &h, const int v_slice,
				const real_1d &angles, const pixel_data &data,
				const real new_h, const real new_angle,
				const int nh, const int na, const int nv);

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

bool CCPi::Nikon_XTek::read_config_file(const std::string path,
					const std::string file)
{
  bool ok = false;
  // read Xtek .xtekct file
  std::string data_file;
  char sep = '\0';
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
	} else if (strncmp(line, "InputSeparator", 14) == 0) {
	  sep = line[15];
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
	  real_1d &h_pixels = set_h_pixels(n_h_pixels);
	  real_1d &v_pixels = set_v_pixels(n_v_pixels);
	  real pixel_base = -((n_h_pixels - 1) * xpixel_size / real(2.0));
	  for (int i = 0; i < n_h_pixels; i++)
	    h_pixels[i] = pixel_base + real(i) * xpixel_size;
	  pixel_base = -((n_v_pixels - 1) * ypixel_size / real(2.0));
	  for (int i = 0; i < n_v_pixels; i++)
	    v_pixels[i] = pixel_base + real(i) * ypixel_size;
	  std::string ctfile;
	  combine_path_and_name(path, "_ctdata.txt", ctfile);
	  ok = read_angles(ctfile, init_angle, n_phi);
	} else {
	  ok = false;
	  report_error("Negative pixel values");
	}
      } else {
	ok = false;
	report_error("Failed to locate all geometry data in xtekct file");
      }
    } else
      report_error("Incorrect header on xtekct file");
    input.close();
  } else
    report_error("Open ", file, " failed");
  if (ok) {
    if (sep != '\0')
      basename += sep;
  }
  return ok;
}

bool CCPi::Nikon_XTek::read_angles(const std::string datafile,
				   const real init_angle, const int n)
{
  // since we read actual angles the answer is probably no
  //std::cout << "Should we use initial angle?\n";
  std::ifstream data(datafile.c_str());
  if (data.good()) {
    real_1d &p = set_phi(n);
    // get first 3 string labels
    std::string t;
    data >> t;
    data >> t;
    data >> t;
    // 3 numbers
    real tmp;
    data >> tmp;
    data >> tmp;
    data >> tmp;
    // now 3 labels for the info we need
    data >> t;
    data >> t;
    data >> t;
    // sets of 3 items terminated by ^M, second is the angle.
    for (int i = 0; i < n; i++) {
      data >> tmp;
      data >> tmp;
      // everything else is radians
      p[i] = real(M_PI) * tmp / real(180.0);
      data >> tmp;
    }
    data.close();
    return true;
  } else {
    report_error("Error opening ctdata file");
    return false;
  }
}

bool CCPi::Nikon_XTek::finish_voxel_geometry(real voxel_origin[3],
					     real voxel_size[3], const int nx,
					     const int ny, const int nz) const
{
  //const voxel_data::size_type *s = voxels.shape();
  real size = real(2.0) * mask_radius / real(nx);
  voxel_size[0] = size;
  voxel_size[1] = size;
  voxel_size[2] = size;
  voxel_origin[0] = -voxel_size[0] * real(nx) / real(2.0) + offset[0];
  voxel_origin[1] = -voxel_size[1] * real(ny) / real(2.0) + offset[1];
  voxel_origin[2] = -voxel_size[2] * real(nz) / real(2.0) + offset[2];
  return true;
}

bool CCPi::Nikon_XTek::read_data_size(const std::string path,
				      const bool phantom)
{
  // phantom already done by setup
  if (phantom) {
    const real_1d &h_pixels = get_h_pixels();
    int n_h_pixels = get_num_h_pixels();
    mask_radius = -get_source_x()
      * std::sin(std::atan(h_pixels[n_h_pixels - 1] /
			   (get_detector_x() - get_source_x())));
  } // else already done by read_config_file
  return true;
}

bool CCPi::Nikon_XTek::read_scans(const std::string path, const int offset,
				  const int block_size, const bool first,
				  const bool phantom)
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

  //real *h_pixels = get_h_pixels();
  //int n_h_pixels = get_num_h_pixels();
  //mask_radius = -get_source_x()
  //* std::sin(std::atan(h_pixels[n_h_pixels - 1] /
  //		 (get_detector_x() - get_source_x())));
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
  voxel_data x(boost::extents[nx][ny][nz], boost::c_storage_order());
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

  pixel_data &pixels = create_pixel_data();
  // perform projection step
  forward_project(pixels, x, image_offset, voxel_size, nx, ny, nz);
  //delete [] x;
  // Todo - could do with adding noise.
  offset[0] = 0.0;
  offset[1] = 0.0;
  offset[2] = 0.0;
  return true;
}

bool CCPi::Nikon_XTek::read_images(const std::string path)
{
  bool ok = true;
  pixel_data &pixels = create_pixel_data();
  //sl_int n_rays = get_data_size();
  std::string pathbase;
  combine_path_and_name(path, basename, pathbase);
  char index[8];
  initialise_progress(get_num_angles(), "Loading data...");
  for (sl_int i = 0; (i < get_num_angles() and ok); i++) {
    snprintf(index, 8, "%04d", int(i + 1));
    std::string name = pathbase + index + ".tif";
    ok = read_tiff(name, pixels, i, get_num_h_pixels(), get_num_v_pixels());
    update_progress(i + 1);
  }
  if (ok) {
    /*
    real max_v = 0.0;
    for (sl_int j = 0; j < n_rays; j++)
      if (max_v < pixel_data[j])
	max_v = pixel_data[j];
    if (max_v > white_level + 0.01)
      std::cout << "Values exceed white level\n";
    else
      max_v = white_level;
    */
    real max_v = real(65535.0);
    // scale and take -ve log, due to exponential extinction in sample.
    for (int i = 0; i < get_num_angles(); i++) {
      for (int j = 0; j < get_num_v_pixels(); j++) {
	for (int k = 0; k < get_num_h_pixels(); k++) {
	  if (pixels[i][j][k] < real(1.0))
	    pixels[i][j][k] = - std::log(real(0.00001) / max_v);
	  else
	    pixels[i][j][k] = - std::log(pixels[i][j][k] / max_v);
	}
      }
    }
    find_centre(get_num_v_pixels() / 2 + 1);
  }
  return ok;
}

pixel_type linear(const real x1, const real x2, const pixel_type f1,
		  const pixel_type f2, const real x)
{
  return (f1 * (x2 - x) + f2  * (x - x1)) / (x2 - x1);
}

pixel_type bilinear(const real x1, const real x2, const real y1, const real y2,
		    const pixel_type f11, const pixel_type f12,
		    const pixel_type f21, const pixel_type f22,
		    const real x, const real y)
{
  return (f11 * (x2 - x) * (y2 - y) + f21 * (x - x1) * (y2 - y)
	  + f12 * (x2 - x) * (y - y1) + f22 * (x - x1) * (y - y1)) /
    ((x2 - x1) * (y2 - y1));
}

pixel_type angles_linear(const int ph1, const real_1d &h, const int v_slice,
			 const real_1d &angles, const pixel_data &data,
			 const real new_h, const real new_angle,
			 const int nh, const int na, const int nv)
{
  int pa1 = 0;
  for (; pa1 < na - 1; pa1++) {
    if (angles[pa1 + 1] > new_angle)
      break;
  }
  if (angles[pa1] == new_angle) {
    return data[pa1][v_slice][ph1];
  } else {
    // interp pa1 to pa1 + 1
    if (pa1 == na - 1)
      return linear(angles[pa1], angles[0] + real(2.0 * M_PI),
		    data[pa1][v_slice][ph1], data[0][v_slice][ph1],
		    new_angle);
    else
      return linear(angles[pa1], angles[pa1 + 1],
		    data[pa1][v_slice][ph1], data[pa1 + 1][v_slice][ph1],
		    new_angle);
  }
}

pixel_type angles_bilinear(const int ph1, const real_1d &h, const int v_slice,
			   const real_1d &angles, const pixel_data &data,
			   const real new_h, const real new_angle,
			   const int nh, const int na, const int nv)
{
  int pa1 = 0;
  for (; pa1 < na - 1; pa1++) {
    if (angles[pa1 + 1] > new_angle)
      break;
  }
  if (angles[pa1] == new_angle) {
    return linear(h[ph1], h[ph1 + 1],
		  data[pa1][ph1][v_slice], data[pa1][ph1 + 1][v_slice], new_h);
  } else {
    // interp pa1 to pa1 + 1
    if (pa1 == na - 1)
      return bilinear(h[ph1], h[ph1 + 1], angles[pa1],
		      angles[0] + real(2.0 * M_PI),
		      data[pa1][ph1][v_slice], data[0][ph1][v_slice],
		      data[pa1][ph1 + 1][v_slice], data[0][ph1 + 1][v_slice],
		      new_h, new_angle);
    else
      return bilinear(h[ph1], h[ph1 + 1], angles[pa1], angles[pa1 + 1],
		      data[pa1][ph1][v_slice], data[pa1 + 1][ph1][v_slice],
		      data[pa1][ph1 + 1][v_slice],
		      data[pa1 + 1][ph1 + 1][v_slice],
		      new_h, new_angle);
  }
}

pixel_type interpolate2D(const real_1d &h, const int v_slice,
			 const real_1d &angles, const pixel_data &data,
			 const real new_h, const real new_angle, const int nh,
			 const int na, const int nv)
{
  int ph1 = 0;
  for (; ph1 < nh - 1; ph1++) {
    if (h[ph1 + 1] > new_h)
      break;
  }
  if (h[ph1] == new_h) {
    // no need to interpolate h - linear interp angles
    return angles_linear(ph1, h, v_slice, angles, data, new_h, new_angle,
			 nh, na, nv);
  } else {
    // interp ph1 to ph1 + 1
    return angles_bilinear(ph1, h, v_slice, angles, data, new_h, new_angle,
			   nh, na, nv);
  }
}

void CCPi::Nikon_XTek::find_centre(const int v_slice)
{
  // converted from matlab find_centre which is based on the method described
  // in T. Liu - "Direct central ray determination in computed microtomography"
  // Optical Engineering, April 2009, 046501
  int n_precs = 5;
  real precision[5] = { 1, 0.1, 0.01, 0.001, 0.0001 };
  real midpoint = 0.0;
  int nv = get_num_v_pixels();
  int nh = get_num_h_pixels();
  const real_1d &h_pixels = get_h_pixels();
  int na = get_num_angles();
  const real_1d &ph = get_phi();
  pixel_data &px = get_pixel_data();
  real distance = get_detector_x() - get_source_x();
  std::vector<real> gamma_i(nh);
  std::vector<real> beta(nh);
  std::vector<real> s2(nh);
  for (int i = 0; i < n_precs; i++) {
    real scor = midpoint - real(10.0) * precision[i];
    real M[21];
    for (int j = 0; j < 21; j++) {
      // angle of each ray relative to theoretical central ray
      // common term arctan(s1/h0) in equations (1) and (2)
      for (int k = 0; k < nh; k++)
	gamma_i[k] = std::atan(h_pixels[k] / distance);
      // angle of assumed centre of rotation to central ray
      // common term arctan(c0/h0) in equations (1) and (2)
      real gamma_c = std::atan((midpoint + (j - 10) * precision[i]) / distance);
      // eqn (1) - Matlab differs from paper in + not -, why?
      for (int k = 0; k < nh; k++)
	beta[k] = real(M_PI) + real(2.0) * (gamma_i[k] - gamma_c);
      // eqn (2)
      for (int k = 0; k < nh; k++)
        s2[k] = distance * std::tan(real(2.0) * gamma_c - gamma_i[k]);
      // do summation on M in eqn (4), only include horizontal values within
      // the data set so normalise each by the count
      int count = 0;
      pixel_type sum = 0;
      for (int a = 0; a < na; a++) {
	for (int k = 0; k < nh; k++) {
	  if (s2[k] >= h_pixels[0] and s2[k] <= h_pixels[nh - 1]) {
	    real alpha_beta = ph[a] + beta[k];
	    // can it be < -2pi or > 4pi?
	    if (alpha_beta < real(0.0))
	      alpha_beta += real(2.0 * M_PI);
	    else if (alpha_beta > real(2.0 * M_PI))
	      alpha_beta -= real(2.0 * M_PI);
	    pixel_type p = interpolate2D(h_pixels, v_slice, ph, px,
					 s2[k], alpha_beta, nh, na, nv);
	    // if (p > 0.0) { ?
	    pixel_type t = px[a][k][v_slice] - p;
	    sum += t * t;
	    count++;
	    //} ?
	  }
	}
      }
      M[j] = sum / (pixel_type)count;
    }
    // minimum value and index
    int ind_m = -1;
    real min_m = 1e20;
    for (int j = 0; j < 21; j++) {
      if (min_m > M[j]) {
	ind_m = j;
	min_m = M[j];
      }
    }
    if (ind_m < 0)
      report_error("Failure in XTek find centre");
    real scor_m = scor + ind_m * precision[i];
    add_output("Precision ");
    add_output(precision[i]);
    add_output(": COR = ");
    add_output(scor_m * get_source_x() / distance);
    add_output(", M = ");
    add_output(min_m);
    send_output();
    midpoint = scor_m;
  }
  // transform centre to required value
  real y_centre = midpoint * get_source_x() / distance;
  set_source(get_source_x(), get_source_y() + y_centre, get_source_z());
  adjust_h_pixels(y_centre);
}

void CCPi::Nikon_XTek::apply_beam_hardening()
{
  // Todo - does this belong in the base class?
  pixel_data &pixels = get_pixel_data();
  for (sl_int i = 0; i < get_num_angles(); i++)
    for (sl_int j = 0; j < get_num_v_pixels(); j++)
      for (sl_int k = 0; k < get_num_h_pixels(); k++)
	pixels[i][j][k] = pixels[i][j][k] * pixels[i][j][k];
}
