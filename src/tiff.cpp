
#include <iostream>
#include "tiffio.h"

#include "base_types.hpp"
#include "tiff.hpp"
#include "timer.hpp"

#ifndef USE_TIMER
#  define USE_TIMER false
#endif // USE_TIMER

bool CCPi::read_tiff(const std::string filename, pixel_type pixel_data[],
		     const int n_h_pixels, const int n_v_pixels)
{
  bool ok = true;
  TIFF *tif = TIFFOpen(filename.c_str(), "r");
  if (tif == 0) {
    ok = false;
    std::cerr << "Error opening " << filename << '\n';
  } else {
    timer ldtime(USE_TIMER);
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
	    } else if ((int)w != n_h_pixels) {
	      std::cerr << "Image width mismatch\n";
	      ok = false;
	    } else {
	      tsize_t strip_size = TIFFStripSize(tif);
	      int nstrips = TIFFNumberOfStrips(tif);
	      char *buf = (char *)_TIFFmalloc(strip_size * nstrips);
	      int offset = 0;
	      int result;
	      for (int count = 0; (count < nstrips and ok); count++) {
		if ((result = TIFFReadEncodedStrip(tif, count, buf + offset,
						   strip_size)) == -1) {
		  std::cerr << "Read error in " << filename << '\n';
		  ok = false;
		} else
		  offset += result;
	      }
	      int extra = offset % 2;
	      offset /= 2;
	      if (extra == 1) {
		std::cerr << "Odd number of bytes for 16 bit image\n";
		ok = false;
	      } else if (offset != n_h_pixels * n_v_pixels) {
		std::cerr << "Image size mismatch\n";
		ok = false;
	      } else {
		// (0,0) at top left means vertical order needs to be reversed
		// for voxel box where the origin is at the bottom
		// It looks like strips are [h][v] col-major
		// storage order is [vert][horiz][angles] fortran order
		// for compatibility with Matlab.
		uint16 *b = (uint16 *)buf;
		long angle_offset = 0;
		for (int h = 0; h < n_h_pixels; h++) {
		  for (int v = 0; v < n_v_pixels; v++) {
		    pixel_data[angle_offset] =
		      pixel_type(b[(n_v_pixels - v - 1) * n_h_pixels + h]);
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
    ldtime.accumulate();
    ldtime.output("Tiff load");
  }
  return ok;
}

bool CCPi::write_tiff(const std::string filename, unsigned short sdata[],
		      unsigned char cdata[], const int nx, const int ny,
		      const int width)
{
  bool ok = true;
  TIFF *tif = TIFFOpen(filename.c_str(), "w");
  if (tif == 0) {
    ok = false;
    std::cerr << "Error opening " << filename << '\n';
  } else {
    int strip_size = nx * width / 8; 
    char *buf = (char *)_TIFFmalloc(strip_size);
    uint16 config = PLANARCONFIG_CONTIG;
    if (TIFFSetField(tif, TIFFTAG_PLANARCONFIG, config) == 0) {
      ok = false;
      std::cerr << "Error setting strip image in tiff\n";
    } else {
      uint16 bps = width;
      uint16 spp = 1;
      if (TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, bps) == 0 or
	  TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, spp) == 0) {
	std::cerr << "TIFF error writing bits per sample\n";
	ok = false;
      } else {
	uint16 photo = PHOTOMETRIC_MINISBLACK;
	if (TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, photo) == 0) {
	  std::cerr << "Error setting TIFF photometric info\n";
	  ok = false;
	} else {
	  uint32 w = nx;
	  uint32 l = ny;
	  if (TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, w) == 0 or
	      TIFFSetField(tif, TIFFTAG_IMAGELENGTH, l) == 0) {
	    std::cerr << "Error setting TIFF size\n";
	    ok = false;
	  } else {
	    uint32 rows_strip = 1;
	    if (TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, rows_strip) == 0) {
	      std::cerr << "Error setting TIFF rows per strip\n";
	      ok = false;
	    } else {
	      int result;
	      for (int count = 0; (count < ny and ok); count++) {
		// We keep the order to match the bottom up orientation
		for (int i = 0; i < strip_size; i++)
		  buf[i] = cdata[(ny - count - 1) * strip_size + i];
		if ((result = TIFFWriteEncodedStrip(tif, count, buf,
						    strip_size)) == -1) {
		  std::cerr << "Write error in " << filename << '\n';
		  ok = false;
		}
	      }
	    }
	  }
	}
      }
    }
    _TIFFfree(buf);
    TIFFClose(tif);
  }
  return ok;
}
