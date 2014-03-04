
#ifndef CCPI_MEX_TYPES
#define CCPI_MEX_TYPES

#include <string>
#include <climits>

#if LONG_MAX == 2147483647L
typedef long long sl_int;
#else
typedef long sl_int;
#endif // LONG_MAX

typedef double real;

typedef float voxel_type;
typedef float pixel_type;

typedef double recon_type;

// dummy
typedef double *voxel_data;

#endif // CCPI_MEX_TYPES
