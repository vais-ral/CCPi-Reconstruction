
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

pixel_data &CCPi::instrument::get_pixel_data()
{
  return *pixels;
}

pixel_data &CCPi::instrument::create_pixel_data()
{
  pixels = new pixel_data(boost::extents[n_angles][n_vertical_pixels][n_horizontal_pixels]);
  return *pixels;
}
