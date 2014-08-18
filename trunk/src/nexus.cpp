
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

bool CCPi::read_NeXus(pixel_data &pixels, pixel_2d &i_dark, pixel_2d &f_dark,
		      pixel_2d &i_bright, pixel_2d &f_bright, int &nh_pixels,
		      int &nv_pixels, std::vector<real> &angles, int &nangles,
		      real &hsize, real &vsize, const std::string filename,
		      const bool all_angles, const bool read_data,
		      const int start_idx, const int block_size)
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
  int n_ibright = 0;
  int n_fbright = 0;
  int n_idark = 0;
  int n_fdark = 0;
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
	sl_int n_ang = 0;
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
	if (!angles_ok) {
	  ok = false;
	  std::cerr << "Problem reading angles\n";
	} else {
	  angles.resize(n_ang);
	  input.closeGroup();
	  // now read the data
	  input.openGroup("instrument", "NXinstrument");
	  entries = input.getEntries();
	  std::string data_name;
	  for (data_entry = entries.begin();
	       data_entry != entries.end(); ++data_entry) {
	    // sample is the official name in the NeXus documentation
	    if (data_entry->second == "NXdetector") {
	      //if (data_entry->first == "bright_field" or
	      //  data_entry->first == "dark_field")
	      //std::cerr << "What to do with bright/dark fields?\n";
	      //else if (data_entry->first == "sample" or
	      if (data_entry->first == "detector")
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
	    bool has_size = true;
	    hsize = 1.0;
	    try {
	      input.openData("x_pixel_size");
	    }
	    catch (...) {
	      has_size = false;
	    }
	    if (has_size) {
	      if (info.type != NeXus::CHAR) {
		std::vector<double> tmp(1);
		try {
		  input.getDataCoerce(tmp);
		}
		catch (NeXus::Exception) {
		  has_size = false;
		}
		if (has_size)
		  hsize = tmp[0];
	      }
	      input.closeData();
	    }
	    has_size = true;
	    vsize = 1.0;
	    try {
	      input.openData("y_pixel_size");
	    }
	    catch (...) {
	      has_size = false;
	    }
	    if (has_size) {
	      if (info.type != NeXus::CHAR) {
		std::vector<double> tmp(1);
		try {
		  input.getDataCoerce(tmp);
		}
		catch (NeXus::Exception) {
		  has_size = false;
		}
		if (has_size)
		  vsize = tmp[0];
	      }
	      input.closeData();
	    } else
	      vsize = 1.0;
	    // get keys to data scans
	    input.openData("image_key");
	    std::vector<double> keys;
	    input.getDataCoerce(keys);
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
	      index[1] = info.dims[1] - (start_idx + block_size);
	      index[2] = 0;
	      sizes[0] = 1;
	      sizes[1] = block_size;
	      sizes[2] = (sl_int)info.dims[2];
	      nh_pixels = sizes[2];
	      if (read_data)
                nv_pixels = sizes[1];
              else
                nv_pixels = info.dims[1];
	      //sl_int offset = sizes[1] * sizes[2];
	      uint16_t *ptr = 0;
	      if (read_data)
		ptr = new uint16_t[sizes[1] * sizes[2]];
	      int n_angles = 0;
	      for (sl_int i = 0; i < n_ang; i++) {
		index[0] = i;
		if (read_data)
		  input.getSlab(ptr, index, sizes);
		// Todo - don't want to depend on pixel order here
		// angles,x,y suggests y(v) varies fastest, since given the
		// way the data is taken it doesn't make much sense for
		// angles to vary fastest, so C storage order I guess
		// should we check for an offset attribute?
		if (keys[i] == 0) {
		  // sample
		  if (read_data) {
		    for (sl_int k = 0; k < sizes[1]; k++) {
		      for (sl_int j = 0; j < sizes[2]; j++) {
			pixels[n_angles][j][(block_size - k - 1)] =
			  pixel_type(ptr[k * sizes[2] + j]);
		      }
		    }
		  }
		  angles[n_angles] = real(M_PI) * angle_data[i] / real(180.0);
		  n_angles++;
		} else if (keys[i] == 1) {
		  if (read_data) {
		    // bright
		    // Todo - check angles, no average? ...
		    pixel_type *bptr = i_bright.data();
		    if (i > n_ang / 2) {
		      bptr = f_bright.data();
		      n_fbright++;
		    } else
		      n_ibright++;
		    sl_int j = 0;
		    for (int h = 0; h < sizes[2]; h++) {
		      for (int v = 0; v < sizes[1]; v++) {
			bptr[j] += pixel_type(ptr[v * sizes[2] + h]);
			j++;
		      }
		    }
		  }
		} else if (keys[i] == 2) {
		  if (read_data) {
		    // dark
		    pixel_type *dptr = i_dark.data();
		    if (i > n_ang / 2) {
		      dptr = f_dark.data();
		      n_fdark++;
		    } else
		      n_idark++;
		    sl_int j = 0;
		    for (int h = 0; h < sizes[2]; h++) {
		      for (int v = 0; v < sizes[1]; v++) {
			dptr[j] += pixel_type(ptr[v * sizes[2] + h]);
			j++;
		      }
		    }
		  }
		}
	      }
	      if (read_data)
		delete [] ptr;
	      else
		nangles = n_angles;
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
  if (read_data) {
    ldtime.accumulate();
    ldtime.output("NeXus load");
    // Average bright/dark frames - Todo, something else? outside here?
    for (int i = 0; i < nh_pixels; i++)
      for (int j = 0; j < nv_pixels; j++)
	i_dark[i][j] /= pixel_type(n_idark);
    for (int i = 0; i < nh_pixels; i++)
      for (int j = 0; j < nv_pixels; j++)
	f_dark[i][j] /= pixel_type(n_fdark);
    for (int i = 0; i < nh_pixels; i++)
      for (int j = 0; j < nv_pixels; j++)
	i_bright[i][j] /= pixel_type(n_ibright);
    for (int i = 0; i < nh_pixels; i++)
      for (int j = 0; j < nv_pixels; j++)
	f_bright[i][j] /= pixel_type(n_fbright);
    // destructor closes file
  }
  return ok;
}
