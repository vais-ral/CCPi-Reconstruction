
#include <iostream>
#include "base_types.hpp"
#include "instruments.hpp"

/*
   Todo - should improve structure of pixel_data and not require the rest
   of the code to know the order, although this may be hard to achieve with the
   current Matlab structures.
*/

long CCPi::instrument::get_data_size() const
{
  return n_angles * n_vertical_pixels * n_horizontal_pixels;
}

pixel_type *const CCPi::instrument::get_pixel_data() const
{
  return pixel_data;
}

pixel_type *CCPi::instrument::create_pixel_data()
{
  long n_rays = n_angles * n_vertical_pixels * n_horizontal_pixels;
  pixel_data = new pixel_type[n_rays];
  return pixel_data;
}

void CCPi::instrument::set_phi(real *p, const int n)
{
  phi = p;
  theta = new real[n];
  for (int i = 0; i < n; i++)
    theta[i] = 0.0;
  n_angles = n;
}
