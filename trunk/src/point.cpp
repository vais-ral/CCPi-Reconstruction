
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
		      const sl_int i, const sl_int j, const sl_int k)
{
    
  int /*N_x, N_y, N_z,*/ im_size_x, im_size_y, im_size_z;
    double b_x, b_y, b_z, d_x, d_y, d_z, d_conv;
    double p1_x, p1_y, p1_z, p2_x, p2_y, p2_z;
    
    int x_defined, y_defined, z_defined;
    
    double alpha_x_min, alpha_y_min, alpha_z_min, alpha_x_max, alpha_y_max, 
      alpha_z_max, alpha_min, alpha_max;
    
    //sl_int ray_index;
	
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

bool test_2D(const real start[], const real end[],
	     const real b_x, const real b_y,
	     const real d_x, const real d_y,
	     const int im_size_x, const int im_size_y,
	     recon_type &length1, recon_type &length2)
{
    
  int N_x, N_y;
  real p1_x, p1_y, p2_x, p2_y;
    
  int x_defined, y_defined;
  sl_int i=0,j=0;
    
    real alpha_x_min, alpha_y_min, alpha_x_max, alpha_y_max, 
      alpha_min, alpha_max;
    //recon_type alpha_x_u = 0.0, alpha_y_u = 0.0;
    //recon_type l_ij;
    //int i_min, j_min, i_max, j_max, n_count, i_u, j_u;
    //sl_int i_step, j_step;
    
    //sl_int ray_index;

    //int z_offset = 0;

    p1_x = start[0];
    p1_y = start[1];
    p2_x = end[0];
    p2_y = end[1];

    length1 = 0.0;
    length2 = 0.0;
    /*
    start[2] = source_z;
    end[2] = det_z[curr_ray_z];
    */

	N_x=im_size_x+1;
	N_y=im_size_y+1;

    /* moving back to original d_conv*/
    //d_conv=std::sqrt( (p1_x-p2_x)*(p1_x-p2_x) + (p1_y-p2_y)*(p1_y-p2_y) + (p1_z-p2_z)*(p1_z-p2_z));


    x_defined =  !(equal_to_precision(p1_x,p2_x,PRECISION));
    y_defined =  !(equal_to_precision(p1_y,p2_y,PRECISION));
	
    if( !x_defined && !y_defined)
	return false;

    if (x_defined) {
	alpha_x_min=std::min(alpha_fn(0, p1_x, p2_x, b_x, d_x), alpha_fn(N_x-1, p1_x, p2_x, b_x, d_x));
	alpha_x_max=std::max(alpha_fn(0, p1_x, p2_x, b_x, d_x), alpha_fn(N_x-1, p1_x, p2_x, b_x, d_x));
    }
    else {
	alpha_x_min=-2;
	alpha_x_max=2;
	i=(int) floor_j( phi(real(0.0), p1_x, p2_x, b_x, d_x));
	if ( i < 0 || i >= im_size_x)
	    return false;
	//alpha_x=2;
	//i_min = 1;
	//i_max = 0;
    }

    if(y_defined) {
	alpha_y_min=std::min(alpha_fn(0, p1_y, p2_y, b_y, d_y), alpha_fn(N_y-1, p1_y, p2_y, b_y, d_y));
	alpha_y_max=std::max(alpha_fn(0, p1_y, p2_y, b_y, d_y), alpha_fn(N_y-1, p1_y, p2_y, b_y, d_y));
    }
    else {
	alpha_y_min=-2;
	alpha_y_max=2;
	j=(int) floor_j( phi(real(0.0), p1_y, p2_y, b_y, d_y));
	if ( j < 0 || j >= im_size_y)
	    return false;
	//alpha_y=2;
	//j_min = 1;
	//j_max = 0;
    }

    		
    alpha_min=std::max(real(0.0),
		       std::max(alpha_x_min, alpha_y_min));
    alpha_max=std::min(real(1.0),
		       std::min(alpha_x_max, alpha_y_max));

    /* if ray intersects voxel grid */
    if (alpha_min < alpha_max) {

      length1 = alpha_max;
      length2 = alpha_min;
      return true;
      /*
      //long nyz = long(im_size_y) * long(im_size_z);
	if (x_defined && p1_x < p2_x) {
	    if (equal_to_precision(alpha_min,alpha_x_min,PRECISION)==1)
		i_min=1;
	    else
		i_min = (int) ceil_j(phi(alpha_min, p1_x, p2_x, b_x, d_x));

	    if (equal_to_precision(alpha_max,alpha_x_max,PRECISION)==1)
		i_max = N_x - 1;
	    else
		i_max = (int) floor_j( phi(alpha_max, p1_x, p2_x, b_x, d_x));

	    alpha_x=alpha_fn(i_min, p1_x, p2_x, b_x, d_x);
	}

	else if (x_defined) {
	    if (equal_to_precision(alpha_min,alpha_x_min,PRECISION)==1)
		i_max=N_x-2;
	    else
		i_max = (int) floor_j(phi(alpha_min, p1_x, p2_x, b_x, d_x));

	    if (equal_to_precision(alpha_max,alpha_x_max,PRECISION)==1)
		i_min = 0;
	    else
		i_min = (int) ceil_j( phi(alpha_max, p1_x, p2_x, b_x, d_x));

	    alpha_x=alpha_fn(i_max, p1_x, p2_x, b_x, d_x);
	}


	if (y_defined && p1_y < p2_y) {
	    if (equal_to_precision(alpha_min,alpha_y_min,PRECISION)==1)
		j_min=1;
	    else
		j_min = (int) ceil_j(phi(alpha_min, p1_y, p2_y, b_y, d_y));


	    if (equal_to_precision(alpha_max, alpha_y_max,PRECISION)==1)
		j_max = N_y - 1;
	    else
		j_max = (int) floor_j( phi(alpha_max, p1_y, p2_y, b_y, d_y));

	    alpha_y=alpha_fn(j_min, p1_y, p2_y, b_y, d_y);
	}

	else if (y_defined) {

	    if (equal_to_precision(alpha_min,alpha_y_min,PRECISION)==1)
		j_max=N_y-2;
	    else
		j_max = (int) floor_j(phi(alpha_min, p1_y, p2_y, b_y, d_y));


	    if (equal_to_precision(alpha_max, alpha_y_max, PRECISION)==1)
		j_min = 0;
	    else
		j_min = (int) ceil_j( phi(alpha_max, p1_y, p2_y, b_y, d_y));

	    alpha_y=alpha_fn(j_max, p1_y, p2_y, b_y, d_y);
	}



	N_p=(i_max - i_min +1) + (j_max - j_min + 1);

	if (x_defined) {
	  i=(int) floor_j( phi( (std::min(alpha_x, alpha_y) + alpha_min)/2, p1_x, p2_x, b_x, d_x) );
	alpha_x_u = d_x/std::abs(p2_x-p1_x);
	  if (i < 0)
	    i = 0;
	  else if (i >= im_size_x)
	    i = im_size_x - 1;
	}

	if (y_defined) {
	  j=(int) floor_j( phi( (std::min(alpha_x, alpha_y) + alpha_min)/2, p1_y, p2_y, b_y, d_y) );
	alpha_y_u = d_y/std::abs(p2_y-p1_y);
	  if (j < 0)
	    j = 0;
	  else if (j >= im_size_y)
	    j = im_size_y - 1;
	}

	if (p1_x < p2_x)
	    i_u=1;
	else
	    i_u=-1;

	if (p1_y < p2_y)
	    j_u=1;
	else
	    j_u=-1;

	alpha_c=alpha_min;
	//ray_index = long(i)*nyz + long(j)*long(im_size_z);
	//i_step = long(i_u) * nyz;
	//j_step = long(j_u) * long(im_size_z);
	//k_step = k_u;


	for (n_count=1; n_count<N_p+1;n_count++) {


	    if (x_defined && alpha_x <= alpha_y) {

	      length = alpha_x - alpha_c;
	      return true;

		if( y_defined && alpha_x == alpha_y) {
		    j += j_u;
		    ray_index += j_step;
		    n_count++;
		    alpha_y += alpha_y_u;
		}

		i += i_u;
		ray_index += i_step;
		alpha_c=alpha_x;
		alpha_x += alpha_x_u;
	    }

	    else if (y_defined) {

	      length = alpha_y - alpha_c;
	      return true;
		j=j+j_u;
		ray_index += j_step;
		alpha_c=alpha_y;
		alpha_y += alpha_y_u;
	    }

	    if( i < 0 || j < 0 || i >= im_size_x || j >= im_size_y)
	      break;
	    


	}
                
	if( (alpha_max - alpha_c) > PRECISION) {
	  // Trap for issues with large umber of threads
	  if( i < 0 || j < 0 || i >= im_size_x || j >= im_size_y) {
#ifdef DEBUG
	    if (alpha_max - alpha_c > recon_type(0.005))
	      std::cerr << "Data left on voxel boundary " << alpha_max - alpha_c << '\n';
#endif //DEBUG
	  } else {
	    l_ij=(alpha_max-alpha_c);
	    
	    length = l_ij;
	    return true;
	  }
	}
	
      */
    } /* of alpha_min < alpha_max */
    
    return false;
}

bool test_3D(const real start[], const real end[],
		const real b_x, const real b_y, const real b_z,
		const real d_x, const real d_y, const real d_z,
		const int im_size_x, const int im_size_y,
		const int im_size_z, recon_type &length)
{
    
  int N_x, N_y, N_z;
    real p1_x, p1_y, p1_z, p2_x, p2_y, p2_z;
    recon_type d_conv;
    
    int x_defined, y_defined, z_defined;
    sl_int i=0,j=0,k=0;
    
    recon_type alpha_x_min, alpha_y_min, alpha_z_min, alpha_x_max, alpha_y_max, 
      alpha_z_max, alpha_min, alpha_max;
    //recon_type alpha_x_u = 0.0, alpha_y_u = 0.0, alpha_z_u = 0.0;
    //recon_type l_ij;
    //int i_min, j_min, k_min, i_max, j_max, k_max, n_count, i_u, j_u, k_u;
    //sl_int i_step, j_step, k_step;
    
    //sl_int ray_index;
	
    p1_x = start[0];
    p1_y = start[1];
    p1_z = start[2];
    p2_x = end[0];
    p2_y = end[1];
    p2_z = end[2];

    //im_size_x = options->im_size_x;
    ///im_size_y = options->im_size_y;
    //im_size_z = options->im_size_z;

	N_x=im_size_x+1;
	N_y=im_size_y+1;
	N_z=im_size_z+1;

	/* d: voxel size */
	//d_x = options->d_x;
	//d_y = options->d_y;
	//d_z = options->d_z;

	/* b: grid offset from origin */
	//b_x = options->b_x;
	//b_y = options->b_y;
	//b_z = options->b_z;



    /* use total lengh=alpha_max-alpha_min instead, to get everage, not sum. */
    /* moving back to original d_conv*/
    d_conv=std::sqrt( (p1_x-p2_x)*(p1_x-p2_x) + (p1_y-p2_y)*(p1_y-p2_y) + (p1_z-p2_z)*(p1_z-p2_z));


    x_defined =  !(equal_to_precision(p1_x,p2_x,PRECISION));
    y_defined =  !(equal_to_precision(p1_y,p2_y,PRECISION));
    z_defined =  !(equal_to_precision(p1_z,p2_z,PRECISION));
	
    if( !x_defined && !y_defined && !z_defined)
	return false;

    if (x_defined) {
	alpha_x_min=std::min(alpha_fn(0, p1_x, p2_x, b_x, d_x), alpha_fn(N_x-1, p1_x, p2_x, b_x, d_x));
	alpha_x_max=std::max(alpha_fn(0, p1_x, p2_x, b_x, d_x), alpha_fn(N_x-1, p1_x, p2_x, b_x, d_x));
    }
    else {
	alpha_x_min=-2;
	alpha_x_max=2;
	i=(int) floor_j( phi(real(0.0), p1_x, p2_x, b_x, d_x));
	if ( i < 0 || i >= im_size_x)
	    return false;
	//alpha_x=2;
	//i_min = 1;
	//i_max = 0;
    }

    if(y_defined) {
	alpha_y_min=std::min(alpha_fn(0, p1_y, p2_y, b_y, d_y), alpha_fn(N_y-1, p1_y, p2_y, b_y, d_y));
	alpha_y_max=std::max(alpha_fn(0, p1_y, p2_y, b_y, d_y), alpha_fn(N_y-1, p1_y, p2_y, b_y, d_y));
    }
    else {
	alpha_y_min=-2;
	alpha_y_max=2;
	j=(int) floor_j( phi(real(0.0), p1_y, p2_y, b_y, d_y));
	if ( j < 0 || j >= im_size_y)
	    return false;
	//alpha_y=2;
	//j_min = 1;
	//j_max = 0;
    }

    		
    if(z_defined) {
	alpha_z_min=std::min(alpha_fn(0, p1_z, p2_z, b_z, d_z), alpha_fn(N_z-1, p1_z, p2_z, b_z, d_z));
	alpha_z_max=std::max(alpha_fn(0, p1_z, p2_z, b_z, d_z), alpha_fn(N_z-1, p1_z, p2_z, b_z, d_z));
    }
    else {
	alpha_z_min=-2;
	alpha_z_max=2;
	k=(int) floor_j( phi(real(0.0), p1_z, p2_z, b_z, d_z));
	if ( k < 0 || k >= im_size_z)
	    return false;
	//alpha_z=2;
	//k_min = 1;
	//k_max = 0;
    }
		
    alpha_min=std::max(real(0.0),
		       max3_dbl(alpha_x_min, alpha_y_min, alpha_z_min));
    alpha_max=std::min(real(1.0),
		       min3_dbl(alpha_x_max, alpha_y_max, alpha_z_max));

    /* if ray intersects voxel grid */
    if (alpha_min < alpha_max) {

      // since its a single voxel we are testing, min/max is the length
      length = (alpha_max - alpha_min) * d_conv;
      return true;
	/*
	if (x_defined && p1_x < p2_x) {
	    if (equal_to_precision(alpha_min,alpha_x_min,PRECISION)==1)
		i_min=1;
	    else
		i_min = (int) ceil_j(phi(alpha_min, p1_x, p2_x, b_x, d_x));

	    if (equal_to_precision(alpha_max,alpha_x_max,PRECISION)==1)
		i_max = N_x - 1;
	    else
		i_max = (int) floor_j( phi(alpha_max, p1_x, p2_x, b_x, d_x));

	    alpha_x=alpha_fn(i_min, p1_x, p2_x, b_x, d_x);
	}

	else if (x_defined) {
	    if (equal_to_precision(alpha_min,alpha_x_min,PRECISION)==1)
		i_max=N_x-2;
	    else
		i_max = (int) floor_j(phi(alpha_min, p1_x, p2_x, b_x, d_x));

	    if (equal_to_precision(alpha_max,alpha_x_max,PRECISION)==1)
		i_min = 0;
	    else
		i_min = (int) ceil_j( phi(alpha_max, p1_x, p2_x, b_x, d_x));

	    alpha_x=alpha_fn(i_max, p1_x, p2_x, b_x, d_x);
	}


	if (y_defined && p1_y < p2_y) {
	    if (equal_to_precision(alpha_min,alpha_y_min,PRECISION)==1)
		j_min=1;
	    else
		j_min = (int) ceil_j(phi(alpha_min, p1_y, p2_y, b_y, d_y));


	    if (equal_to_precision(alpha_max, alpha_y_max,PRECISION)==1)
		j_max = N_y - 1;
	    else
		j_max = (int) floor_j( phi(alpha_max, p1_y, p2_y, b_y, d_y));

	    alpha_y=alpha_fn(j_min, p1_y, p2_y, b_y, d_y);
	}

	else if (y_defined) {

	    if (equal_to_precision(alpha_min,alpha_y_min,PRECISION)==1)
		j_max=N_y-2;
	    else
		j_max = (int) floor_j(phi(alpha_min, p1_y, p2_y, b_y, d_y));


	    if (equal_to_precision(alpha_max, alpha_y_max, PRECISION)==1)
		j_min = 0;
	    else
		j_min = (int) ceil_j( phi(alpha_max, p1_y, p2_y, b_y, d_y));

	    alpha_y=alpha_fn(j_max, p1_y, p2_y, b_y, d_y);
	}



	if (z_defined && p1_z < p2_z) {
	    if (equal_to_precision(alpha_min,alpha_z_min,PRECISION)==1)
		k_min=1;
	    else
		k_min = (int) ceil_j(phi(alpha_min, p1_z, p2_z, b_z, d_z));


	    if (equal_to_precision(alpha_max, alpha_z_max,PRECISION)==1)
		k_max = N_z - 1;
	    else
		k_max = (int) floor_j( phi(alpha_max, p1_z, p2_z, b_z, d_z));

	    alpha_z=alpha_fn(k_min, p1_z, p2_z, b_z, d_z);
	}

	else if (z_defined) {

	    if (equal_to_precision(alpha_min,alpha_z_min,PRECISION)==1)
		k_max=N_z-2;
	    else
		k_max = (int) floor_j(phi(alpha_min, p1_z, p2_z, b_z, d_z));


	    if (equal_to_precision(alpha_max, alpha_z_max, PRECISION)==1)
		k_min = 0;
	    else
		k_min = (int) ceil_j( phi(alpha_max, p1_z, p2_z, b_z, d_z));

	    alpha_z=alpha_fn(k_max, p1_z, p2_z, b_z, d_z);
	}


	N_p=(i_max - i_min +1) + (j_max - j_min + 1) + (k_max - k_min + 1);

	if (x_defined) {
	    i=(int) floor_j( phi( (min3_dbl(alpha_x, alpha_y, alpha_z) + alpha_min)/2, p1_x, p2_x, b_x, d_x) );
	alpha_x_u = d_x/std::abs(p2_x-p1_x);
	  if (i < 0)
	    i = 0;
	  else if (i >= im_size_x)
	    i = im_size_x - 1;
	}

	if (y_defined) {
	    j=(int) floor_j( phi( (min3_dbl(alpha_x, alpha_y, alpha_z) + alpha_min)/2, p1_y, p2_y, b_y, d_y) );
	alpha_y_u = d_y/std::abs(p2_y-p1_y);
	  if (j < 0)
	    j = 0;
	  else if (j >= im_size_y)
	    j = im_size_y - 1;
	}
	if (z_defined) {
	    k=(int) floor_j( phi( (min3_dbl(alpha_x, alpha_y, alpha_z) + alpha_min)/2, p1_z, p2_z, b_z, d_z) );
	alpha_z_u = d_z/std::abs(p2_z-p1_z);
	  if (k < 0)
	    k = 0;
	  else if (k >= im_size_z)
	    k = im_size_z - 1;
	}

	if (p1_x < p2_x)
	    i_u=1;
	else
	    i_u=-1;

	if (p1_y < p2_y)
	    j_u=1;
	else
	    j_u=-1;

	if (p1_z < p2_z)
	    k_u=1;
	else
	    k_u=-1;


	alpha_c=alpha_min;
	//ray_index = i*im_size_y*im_size_z + j*im_size_z + (k+z_offset);
	//i_step = i_u * im_size_y * im_size_z;
	//j_step = j_u * im_size_z;
	//k_step = k_u;

	for (n_count=1; n_count<N_p+1;n_count++) {


	    if (x_defined && alpha_x <= alpha_y && alpha_x <= alpha_z) {

	      length = (alpha_x - alpha_c) * d_conv;
	      return true;

		if( y_defined && alpha_x == alpha_y) {
		    j += j_u;
		    //ray_index += j_step;
		    n_count++;
		    alpha_y += alpha_y_u;
		}

		if( z_defined && alpha_x == alpha_z) {
		    k += k_u;
		    //ray_index += k_step;
		    n_count++;
		    alpha_z += alpha_z_u;
		}

		i += i_u;
		//ray_index += i_step;
		alpha_c=alpha_x;
		alpha_x += alpha_x_u;
	    }

	    else if (y_defined && alpha_y <= alpha_z) {

	      length = (alpha_y - alpha_c) * d_conv;
	      return true;

		if( z_defined && alpha_y == alpha_z) {
		    k += k_u;
		    //ray_index += k_step;
		    n_count++;
		    alpha_z += alpha_z_u;
		}

		j=j+j_u;
		//ray_index += j_step;
		alpha_c=alpha_y;
		alpha_y += alpha_y_u;
	    }

	    else if (z_defined) {

	      length = (alpha_z - alpha_c) * d_conv;
	      return true;

		k += k_u;
		//ray_index += k_step;
		alpha_c=alpha_z;
		alpha_z += alpha_z_u;
	    }


	    if( i < 0 || j < 0 || k < 0 || i >= im_size_x || j >= im_size_y || k >= im_size_z)
	      break;
	    


	}
                
	if( (alpha_max - alpha_c) > PRECISION) {
	  // Trap for issues with large umber of threads
	  if( i < 0 || j < 0 || k < 0 || i >= im_size_x || j >= im_size_y || k >= im_size_z) {
#ifdef DEBUG
	    if (alpha_max - alpha_c > recon_type(0.005))
	      std::cerr << "Data left on voxel boundary " << alpha_max - alpha_c << '\n';
#endif //DEBUG
	  } else {
	    //l_ij=(alpha_max-alpha_c)*d_conv;

	    length = (alpha_max - alpha_c) * d_conv;
	    return true;
	  }
	}
	*/
	
    } /* of alpha_min < alpha_max */
    
    return false;
}
