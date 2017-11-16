
#ifndef MPI_WRAPPER
#define MPI_WRAPPER
#include "CCPiDefines.h"
namespace machine {

  CCPI_EXPORT void __cdecl initialise(const int nthreads = 0);
  CCPI_EXPORT void __cdecl exit();
  CCPI_EXPORT int __cdecl get_number_of_processors();
  CCPI_EXPORT int __cdecl get_processor_id();

}

#endif // MPI_WRAPPER
