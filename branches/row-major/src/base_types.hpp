
#ifndef CCPI_RECON_TYPES
#define CCPI_RECON_TYPES

#include <climits>
#include <string>
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

typedef boost::multi_array<voxel_type, 3> voxel_data;

#endif // CCPI_RECON_TYPES
