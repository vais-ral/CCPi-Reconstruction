
#ifndef CCPI_RECON_TYPES
#define CCPI_RECON_TYPES

#include <string>
#include <boost/multi_array.hpp>

typedef double real;
typedef double voxel_type;
typedef double pixel_type;
typedef double recon_type;

typedef boost::multi_array<voxel_type, 3> voxel_data;

#endif // CCPI_RECON_TYPES
