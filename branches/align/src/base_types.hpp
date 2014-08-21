
#ifndef CCPI_RECON_TYPES
#define CCPI_RECON_TYPES

#include <climits>
#include <string>
#include <vector>
#include <boost/multi_array.hpp>
#include "allocator.hpp"

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

typedef std::vector<real, aligned_allocator<real> > real_1d;
typedef std::vector<real, aligned_allocator<real> > real_1dr;
typedef std::vector<pixel_type, aligned_allocator<pixel_type> > pixel_1d;
typedef std::vector<voxel_type, aligned_allocator<voxel_type> > voxel_1d;

typedef boost::multi_array<real, 2, aligned_allocator<real> > real_2d;
typedef boost::multi_array<real, 3, aligned_allocator<real> > real_3d;
typedef boost::multi_array<pixel_type, 2,
			   aligned_allocator<pixel_type> > pixel_2d;
typedef boost::multi_array<pixel_type, 3,
			   aligned_allocator<pixel_type> > pixel_3d;
typedef boost::multi_array<voxel_type, 2,
			   aligned_allocator<voxel_type> > voxel_2d;
typedef boost::multi_array<voxel_type, 3,
			   aligned_allocator<voxel_type> > voxel_3d;

typedef std::vector<recon_type, aligned_allocator<recon_type> > recon_1d;
typedef boost::multi_array<recon_type, 2,
			   aligned_allocator<recon_type> > recon_2d;

typedef boost::multi_array<pixel_type, 3,
			   aligned_allocator<pixel_type> > pixel_data;
typedef boost::multi_array<voxel_type, 3,
			   aligned_allocator<voxel_type> > voxel_data;

typedef std::vector<int, aligned_allocator<int> > int_1d;
typedef std::vector<long, aligned_allocator<long> > long_1d;

#endif // CCPI_RECON_TYPES
