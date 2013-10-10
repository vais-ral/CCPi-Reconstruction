
#define PRECISION 0.00000001

#include <iostream>
#include <cstdlib>
#include <cmath>
#include "base_types.hpp"

static inline double alpha_fn(const int n, const double p1, const double p2,
			      const double b, const double d)
{
    return ( (b+n*d) - p1)/(p2-p1);
}

static inline double p(const double alpha, const double p1, const double p2)
{
    return p1+alpha*(p2-p1);
}

static inline double phi(const double alpha, const double p1, const double p2,
			 const double b, const double d)
{
    return ( p(alpha, p1, p2)-b)/d;
}

static inline int equal_to_precision(const double x, const double y,
				     const double prec)
{
    return fabs(x-y) < prec;
}

static inline double min_dbl(const double a, const double b)
{
  return a < b ? a : b;
}

static inline double max_dbl(const double a, const double b)
{
  return a > b ? a : b;
}

static inline double min3_dbl(const double a, const double b, const double c)
{
  return a < b ? min_dbl(a,c) : min_dbl(b,c);
}

static inline double max3_dbl(const double a, const double b, const double c)
{
  return a > b ? max_dbl(a,c) : max_dbl(b,c);
}

static inline double ceil_j(const double arg)
{
  return arg == (int)arg ? arg+1 : ceil( arg );
}

static inline double floor_j(const double arg)
{
  return floor( arg );
}

double generate_point(const double start[], const double end[],
		      const int vx, const int vy, const int vz, 
		      const double vsize[3], const double vorigin[3],
		      const long i, const long j, const long k)
{
    
  int /*N_x, N_y, N_z,*/ im_size_x, im_size_y, im_size_z;
    double b_x, b_y, b_z, d_x, d_y, d_z, d_conv;
    double p1_x, p1_y, p1_z, p2_x, p2_y, p2_z;
    
    int x_defined, y_defined, z_defined;
    
    double alpha_x_min, alpha_y_min, alpha_z_min, alpha_x_max, alpha_y_max, 
      alpha_z_max, alpha_min, alpha_max;
    
    //long ray_index;
	
    p1_x = start[0];
    p1_y = start[1];
    p1_z = start[2];
    p2_x = end[0];
    p2_y = end[1];
    p2_z = end[2];

    im_size_x = vx;
    im_size_y = vy;
    im_size_z = vz;

    if (i < 0 || i >= im_size_x || j < 0 || j >= im_size_y ||
	k < 0 || k >= im_size_z) {
      std::cerr << "i,j,k outside voxels\n";
      return -1.0;
    }

    //N_x=im_size_x+1;
    //N_y=im_size_y+1;
    //N_z=im_size_z+1;

    /* d: voxel size */
    d_x = vsize[0];
    d_y = vsize[1];
    d_z = vsize[2];

    /* b: grid offset from origin */
    b_x = vorigin[0];
    b_y = vorigin[1];
    b_z = vorigin[2];

    /* use total lengh=alpha_max-alpha_min instead, to get everage, not sum. */
    /* moving back to original d_conv*/
    d_conv=sqrt( (p1_x-p2_x)*(p1_x-p2_x) + (p1_y-p2_y)*(p1_y-p2_y) + (p1_z-p2_z)*(p1_z-p2_z));


    x_defined =  !(equal_to_precision(p1_x,p2_x,PRECISION));
    y_defined =  !(equal_to_precision(p1_y,p2_y,PRECISION));
    z_defined =  !(equal_to_precision(p1_z,p2_z,PRECISION));
	
    if( !x_defined && !y_defined && !z_defined)
	return 0.0;

    if (x_defined) {
	alpha_x_min=min_dbl(alpha_fn(i, p1_x, p2_x, b_x, d_x), alpha_fn(i+1, p1_x, p2_x, b_x, d_x));
	alpha_x_max=max_dbl(alpha_fn(i, p1_x, p2_x, b_x, d_x), alpha_fn(i+1, p1_x, p2_x, b_x, d_x));
    }
    else {
      if (b_x + i * d_x >= p1_x || b_x + (i + 1) * d_x < p1_x)
	return 0.0;
	alpha_x_min=-2;
	alpha_x_max=2;
    }

    if(y_defined) {
	alpha_y_min=min_dbl(alpha_fn(j, p1_y, p2_y, b_y, d_y), alpha_fn(j+1, p1_y, p2_y, b_y, d_y));
	alpha_y_max=max_dbl(alpha_fn(j, p1_y, p2_y, b_y, d_y), alpha_fn(j+1, p1_y, p2_y, b_y, d_y));
    }
    else {
      if (b_y + j * d_y >= p1_y || b_y + (j + 1) * d_y < p1_y)
	return 0.0;
	alpha_y_min=-2;
	alpha_y_max=2;
    }

    		
    if(z_defined) {
	alpha_z_min=min_dbl(alpha_fn(k, p1_z, p2_z, b_z, d_z), alpha_fn(k+1, p1_z, p2_z, b_z, d_z));
	alpha_z_max=max_dbl(alpha_fn(k, p1_z, p2_z, b_z, d_z), alpha_fn(k+1, p1_z, p2_z, b_z, d_z));
    }
    else {
      if (b_z + k * d_z >= p1_z || b_z + (k + 1) * d_z < p1_z)
	return 0.0;
	alpha_z_min=-2;
	alpha_z_max=2;
    }
		
    alpha_min=max_dbl(0.0, max3_dbl(alpha_x_min, alpha_y_min, alpha_z_min));
    alpha_max=min_dbl(1.0, min3_dbl(alpha_x_max, alpha_y_max, alpha_z_max));

    /* if ray intersects voxel grid */
    if (alpha_min < alpha_max) {

      //ray_index = k*im_size_y*im_size_x + j*im_size_x + i;
        //data = (alpha_max-alpha_min)*d_conv * vol_data[ray_index];

	//*ray_data += (float)data;

      return (alpha_max-alpha_min)*d_conv;
	
    } /* of alpha_min < alpha_max */
    else
      return 0.0;
    
}

void test_voxel(const int i, const int j, const int k,
		const int a, const int v, const int h, const double d,
		const double l,
		const double voxel_origin[3], const double voxel_size[3],
		const double v_pos, const double h_pos, const double cphi,
		const double sphi, const int nx, const int ny, const int nz,
		const char message[])
{
  double det_x = 2.0 * std::max(std::abs(voxel_origin[0]),
				std::max(std::abs(voxel_origin[1]),
					 std::abs(voxel_origin[2])));
  double start[3];
  double end[3];
  end[2] = v_pos;
  start[2] = end[2];
  end[0] = cphi * det_x - sphi * h_pos;
  end[1] = sphi * det_x + cphi * h_pos;
  start[0] = end[0] - 3.0 * cphi * det_x;
  start[1] = end[1] - 3.0 * sphi * det_x;
  double value = generate_point(start, end, nx, ny, nz, voxel_size,
				voxel_origin, i, j, k);
  if (std::abs(d * l - value) >= 1e-5)
    std::cerr << message << ' ' << i << ' ' << j << ' ' << k <<
      ", " << a << ' ' << v << ' ' << h << ", " << d << ' ' << l << ", " <<
      value << ", " << cphi << ' ' << sphi << '\n';
}

void test_voxel(const int i, const int j, const int k,
		const int a, const int v, const int h, const double d,
		const double l,
		const double voxel_origin[3], const double voxel_size[3],
		const double v_pos, const double h_pos1, const double h_pos2,
		const double cphi,
		const double sphi, const int nx, const int ny, const int nz,
		const char message[])
{
  double det_x = 2.0 * std::max(std::abs(voxel_origin[0]),
				std::max(std::abs(voxel_origin[1]),
					 std::abs(voxel_origin[2])));
  double start[3];
  double end[3];
  end[2] = v_pos;
  start[2] = end[2];
  end[0] = cphi * det_x - sphi * h_pos1;
  end[1] = sphi * det_x + cphi * h_pos1;
  start[0] = end[0] - 3.0 * cphi * det_x;
  start[1] = end[1] - 3.0 * sphi * det_x;
  double value = generate_point(start, end, nx, ny, nz, voxel_size,
				voxel_origin, i, j, k);
  end[0] = cphi * det_x - sphi * h_pos2;
  end[1] = sphi * det_x + cphi * h_pos2;
  start[0] = end[0] - 3.0 * cphi * det_x;
  start[1] = end[1] - 3.0 * sphi * det_x;
  value += generate_point(start, end, nx, ny, nz, voxel_size,
			  voxel_origin, i, j, k);
  if (std::abs(d * l - value) >= 1e-5)
    std::cerr << message << ' ' << i << ' ' << j << ' ' << k <<
      ", " << a << ' ' << v << ' ' << h << ", " << d << ' ' << l << ", " <<
      value << ", " << cphi << ' ' << sphi << '\n';
}
