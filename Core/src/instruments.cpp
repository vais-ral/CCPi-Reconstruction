
#include <cmath>
#include <iostream>
#ifdef MATLAB_MEX_FILE
#  include "mex_types.hpp"
#else
#  include "base_types.hpp"
#endif // mex
#include "instruments.hpp"
#include "ui_calls.hpp"

/*
   Todo - should improve structure of pixel_data and not require the rest
   of the code to know the order, although this may be hard to achieve with the
   current Matlab structures.
*/

CCPi::instrument::~instrument()
{
  if (pixels != 0)
    delete pixels;
}

pixel_3d &CCPi::instrument::get_pixel_data()
{
  return *pixels;
}

pixel_3d &CCPi::instrument::create_pixel_data()
{
  pixels = new pixel_3d(boost::extents[n_angles][n_horizontal_pixels][n_vertical_pixels]);
  return *pixels;
}

int CCPi::instrument::calc_v_alignment(const int n, const int pix_per_vox,
				       const bool cone)
{
  data_v_size = n;
  int nvox = n / pix_per_vox;
  if (n % pix_per_vox != 0)
    nvox++;
  int align = aligned_allocator<voxel_type>::alignment / sizeof(voxel_type);
  int vox_blocks = nvox / align;
  if (nvox % align != 0)
    vox_blocks++;
  // make sure the cone can split into 2 equal aligned halves
  if (cone) {
    if (vox_blocks %2 != 0)
      vox_blocks++;
  }
  int npix = vox_blocks * align * pix_per_vox;
  int diff = npix - n;
  data_v_offset = diff / 2;
  // try and check that midpoint is in aligned middle of pixels?
  // if n is an even number then they must be, mustn't they?
  /*
  add_output("Size is ");
  add_output(npix);
  add_output(" ");
  add_output(data_v_offset);
  send_output();
  */
  return npix;
}
