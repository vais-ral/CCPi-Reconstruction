
#ifndef CCPI_MEX_TYPES
#define CCPI_MEX_TYPES

#include <string>
#include <climits>
#include <vector>
#include <boost/multi_array.hpp>

#if LONG_MAX == 2147483647L
typedef long long sl_int;
#else
typedef long sl_int;
#endif // LONG_MAX

typedef double real;

typedef float voxel_type;
typedef float pixel_type;

typedef float recon_type;

typedef std::vector<real> real_1d;
typedef boost::multi_array_ref<real, 1> real_1dr;

typedef std::vector<recon_type> recon_1d;
typedef boost::multi_array<recon_type, 2> recon_2d;

typedef boost::multi_array_ref<pixel_type, 3> pixel_data;
typedef boost::multi_array<pixel_type, 3> pixel_3d;
typedef boost::multi_array_ref<voxel_type, 3> voxel_data;
typedef boost::multi_array<voxel_type, 3> voxel_3d;

#endif // CCPI_MEX_TYPES
