
#ifndef CCPI_RECON_TYPES
#define CCPI_RECON_TYPES

#include <climits>
#include <string>
#include <vector>
#include <boost/multi_array.hpp>

#ifdef WIN32
#  include <ciso646>
#  define snprintf _snprintf
#endif // WIN32

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
typedef std::vector<real> real_1dr;
typedef std::vector<pixel_type> pixel_1d;
typedef std::vector<voxel_type> voxel_1d;

typedef boost::multi_array<real, 2> real_2d;
typedef boost::multi_array<real, 3> real_3d;
typedef boost::multi_array<pixel_type, 2> pixel_2d;
typedef boost::multi_array<pixel_type, 3> pixel_3d;
typedef boost::multi_array<voxel_type, 2> voxel_2d;
typedef boost::multi_array<voxel_type, 3> voxel_3d;

typedef std::vector<recon_type> recon_1d;
typedef boost::multi_array<recon_type, 2> recon_2d;

typedef boost::multi_array<pixel_type, 3> pixel_data;
typedef boost::multi_array<voxel_type, 3> voxel_data;

#endif // CCPI_RECON_TYPES
