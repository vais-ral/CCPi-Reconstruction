
#ifdef MATLAB_MEX_FILE
#  ifdef WIN32
#    define snprintf _snprintf
#  endif // WIN32
#  include "mex_types.hpp"
#else
#  include "base_types.hpp"
#endif // MEX_FILE
#include "mpi.hpp"
#include "omp.h"
#include "accel.hpp"

#ifdef USE_MPI

#  include "mpi.h"

void machine::initialise(const int nthreads)
{
  int argc = 0;
  char **argv = 0;
  MPI_Init(&argc, &argv);
  if (nthreads > 0)
    omp_set_num_threads(nthreads);
  init_accelerator();
}

void machine::exit()
{
  close_accelerator();
  MPI_Finalize();
}

int machine::get_number_of_processors()
{
  int num_nodes;
  MPI_Comm_size(MPI_COMM_WORLD, &num_nodes);
  return num_nodes;
}

int machine::get_processor_id()
{
  int node_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &node_id);
  return node_id;
}

#else

void machine::initialise(const int nthreads)
{
  if (nthreads > 0)
    omp_set_num_threads(nthreads);
  init_accelerator();
}

void machine::exit()
{
  close_accelerator();
}

int machine::get_number_of_processors()
{
  return 1;
}

int machine::get_processor_id()
{
  return 0;
}

#endif // USE_MPI
