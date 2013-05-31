
#include <iostream>
#include <string>
#include <map>
#include "NeXusFile.hpp"
#include "NeXusException.hpp"
#include "base_types.hpp"
#include "nexus.hpp"
#include "timer.hpp"

#ifndef USE_TIMER
#  define USE_TIMER false
#endif // USE_TIMER

bool CCPi::read_NeXus(pixel_type * &pixels, int &nh_pixels, int &nv_pixels,
		      real * &angles, int &nangles, real &hsize, real &vsize,
		      const std::string filename, const bool all_angles)
{
  bool ok = true;
  timer ldtime(USE_TIMER);
  NeXus::File input(filename);
  // catch exception if open/construct failed?
  /*
  bool ok = true;
  try {
    input.openGroup("entry", "NXentry");
  }
  catch (NeXus::Exception &excp) {
    ok = false;
    // Dumb library writes error whether I trap this or not - Yuck.
    std::cerr << excp.what() << '\n';
  }
  if (ok) {
  */
  std::map<std::string, std::string> entries = input.getEntries();
  std::map<std::string, std::string>::const_iterator data_entry;
  // check for entry:NXentry
  for (data_entry = entries.begin();
       data_entry != entries.end(); ++data_entry) {
    if (data_entry->second == "NXentry")
      if (data_entry->first == "entry")
	break;
  }
  if (data_entry == entries.end()) {
    // retry with entry:NXentry*
    for (data_entry = entries.begin();
	 data_entry != entries.end(); ++data_entry) {
      if (data_entry->second == "NXentry")
	if (data_entry->first.compare(0, 5, "entry") == 0)
	  break;
    }
  }
  if (data_entry == entries.end()) {
    ok = false;
    std::cerr << "Couldn't find NeXus file main entry\n";
  } else {
    input.openGroup(data_entry->first, "NXentry");
    entries = input.getEntries();
    // look for the tomographic data in the file
    for (data_entry = entries.begin();
	 data_entry != entries.end(); ++data_entry) {
      if (data_entry->second == "NXtomo")
	break;
    }
    if (data_entry == entries.end()) {
      for (data_entry = entries.begin();
	   data_entry != entries.end(); ++data_entry) {
	if (data_entry->second == "NXsubentry")
	  if (data_entry->first == "tomo_entry")
	    break;
      }
    }
    if (data_entry == entries.end()) {
      ok = false;
      std::cerr << "Can't find tomographic data\n";
    } else {
      input.openGroup(data_entry->first, data_entry->second);
      // check its NXtomo
      NeXus::Info info;
      input.openData("definition");
      if (info.type == NeXus::CHAR and info.dims.size() == 1) {
	if (input.getStrData() != "NXtomo") {
	  ok = false;
	  std::cerr << "Not tomographic data\n";
	}
      }
      input.closeData();
      if (ok) {
	// open the sample - Todo trap errors?
	input.openGroup("sample", "NXsample");
	// check for sample shifts we can't handle yet
	input.openData("x_translation");
	info = input.getInfo();
	if (info.type == NeXus::FLOAT32 or info.type == NeXus::FLOAT64) {
	  if (info.dims.size() == 1) {
	    std::vector<double> result;
	    input.getDataCoerce(result);
	    bool change = false;
	    for (int64_t i = 1; i < info.dims[0]; i++)
	      change = change or (result[i] != result[0]);
	  }
	}
	input.closeData();
	input.openData("y_translation");
	info = input.getInfo();
	if (info.type == NeXus::FLOAT32 or info.type == NeXus::FLOAT64) {
	  if (info.dims.size() == 1) {
	    std::vector<double> result;
	    input.getDataCoerce(result);
	    bool change = false;
	    for (int64_t i = 1; i < info.dims[0]; i++)
	      change = change or (result[i] != result[0]);
	  }
	}
	input.closeData();
	input.openData("z_translation");
	info = input.getInfo();
	if (info.type == NeXus::FLOAT32 or info.type == NeXus::FLOAT64) {
	  if (info.dims.size() == 1) {
	    std::vector<double> result;
	    input.getDataCoerce(result);
	    bool change = false;
	    for (int64_t i = 1; i < info.dims[0]; i++)
	      change = change or (result[i] != result[0]);
	  }
	}
	input.closeData();
	// get the angles
	bool angles_ok = true;
	long n_ang = 0;
	std::vector<double> angle_data;
	input.openData("rotation_angle");
	info = input.getInfo();
	if (info.type == NeXus::FLOAT32 or info.type == NeXus::FLOAT64) {
	  if (info.dims.size() == 1) {
	    n_ang = info.dims[0];
	    input.getDataCoerce(angle_data);
	  } else
	    angles_ok = false;
	} else
	  angles_ok = false;
	input.closeData();
	long min_angle = 0;
	long max_angle = n_ang;
	if (!angles_ok) {
	  ok = false;
	  std::cerr << "Problem reading angles\n";
	} else {
	  if (!all_angles) {
	    for (long i = 1; i < n_ang; i++)
	      if (angle_data[i] == angle_data[min_angle])
		min_angle++;
	    for (long i = n_ang - 1; i > 1; i--)
	      if (angle_data[i - 1] == angle_data[max_angle - 1])
		max_angle--;
	  }
	  /*
	    std::cout << "Got " << n_ang << " angles\n";
	    std::cout << "from " << angle_data[min_angle] << ' '
	    << min_angle << '\n';
	    std::cout << "  to " << angle_data[max_angle - 1] << ' '
	    << max_angle - 1 << '\n';
	    std::cout << angle_data[19] << ' ';
	    std::cout << angle_data[20] << ' ';
	    std::cout << angle_data[21] << '\n';
	    std::cout << angle_data[1819] << ' ';
	    std::cout << angle_data[1820] << ' ';
	    std::cout << angle_data[1821] << '\n';
	  */
	  angles = new real[max_angle - min_angle];
	  for (int i = min_angle; i < max_angle; i++)
	    angles[i - min_angle] = angle_data[i];
	  nangles = max_angle - min_angle;
	  input.closeGroup();
	  // now read the data
	  input.openGroup("instrument", "NXinstrument");
	  entries = input.getEntries();
	  std::string data_name;
	  for (data_entry = entries.begin();
	       data_entry != entries.end(); ++data_entry) {
	    // sample is the official name in the NeXus documentation
	    if (data_entry->second == "NXdetector") {
	      if (data_entry->first == "bright_field" or
		  data_entry->first == "dark_field")
		std::cerr << "What to do with bright/dark fields?\n";
	      else if (data_entry->first == "sample" or
		       data_entry->first == "detector")
		data_name = data_entry->first;
	    }
	  }
	  if (data_name == "") {
	    ok = false;
	    std::cerr << "Couldn't find actual data\n";
	  } else {
	    input.openGroup(data_name, "NXdetector");
	    // Todo - get distance?
	    // pixel sizes - Todo calibrate scale um/mm/m...
	    input.openData("x_pixel_size");
	    std::vector<double> tmp(1);
	    input.getDataCoerce(tmp);
	    hsize = tmp[0];
	    input.closeData();
	    input.openData("y_pixel_size");
	    input.getDataCoerce(tmp);
	    vsize = tmp[0];
	    input.closeData();
	    // read data
	    input.openData("data");
	    info = input.getInfo();
	    if (info.type != NeXus::UINT16) {
	      ok = false;
	      std::cerr << "Unexpected data type\n";
	    } else if (info.dims.size() != 3 or info.dims[0] != n_ang) {
	      ok = false;
	      std::cerr << "Unexpected data range\n";
	    } else {
	      std::vector<long> index(3);
	      std::vector<long> sizes(3);
	      index[0] = 0;
	      index[1] = 0;
	      index[2] = 0;
	      sizes[0] = 1;
	      sizes[1] = (long)info.dims[1];
	      sizes[2] = (long)info.dims[2];
	      nh_pixels = sizes[1];
	      nv_pixels = sizes[2];
	      pixels = new pixel_type[long(nangles) * long(sizes[1])
				      * long(sizes[2])];
	      uint16_t *ptr = new uint16_t[sizes[1] * sizes[2]];
	      long offset = sizes[1] * sizes[2];
	      for (long i = min_angle; i < max_angle; i++) {
		index[0] = i;
		input.getSlab(ptr, index, sizes);
		// Todo - don't want to depend on pixel order here
		// angles,x,y suggests y(v) varies fastest, since given the
		// way the data is taken it doesn't make much sense for
		// angles to vary fastest, so C storage order I guess
		// should we check for an offset attribute?
		for (long j = 0; j < sizes[1]; j++) {
		  for (long k = 0; k < sizes[2]; k++) {
		    pixels[j + k * sizes[1] + (i - min_angle) * offset] =
		      pixel_type(ptr[j * sizes[2] + k]);
		  }
		}
	      }
	      delete [] ptr;
	    }
	    input.closeData();
	    input.closeGroup();
	  }
	}
	input.closeGroup();
      } // !ok
    }
    input.closeGroup();
  }
  ldtime.accumulate();
  ldtime.output("NeXus load");
  // destructor closes file
  return ok;
}
