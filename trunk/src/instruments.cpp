
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

sl_int CCPi::instrument::get_data_size() const
{
  return sl_int(n_angles) * sl_int(n_vertical_pixels)
    * sl_int(n_horizontal_pixels);
}

pixel_type *const CCPi::instrument::get_pixel_data() const
{
  return pixel_data;
}

pixel_type *CCPi::instrument::create_pixel_data()
{
  sl_int n_rays = sl_int(n_angles) * sl_int(n_vertical_pixels)
    * sl_int(n_horizontal_pixels);
  pixel_data = new pixel_type[n_rays];
  return pixel_data;
}

void CCPi::instrument::set_pixel_data(pixel_type *p, const sl_int n)
{
  sl_int n_rays = sl_int(n_angles) * sl_int(n_vertical_pixels)
    * sl_int(n_horizontal_pixels);
  if (n != n_rays)
    report_error("Size mismatch setting pixel data");
  pixel_data = p;
}

/*
void CCPi::instrument::set_phi(real *p, const int n)
{
  phi = p;
  n_angles = n;
}
*/
