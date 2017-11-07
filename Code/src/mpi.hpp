
#ifndef MPI_WRAPPER
#define MPI_WRAPPER

namespace machine {

  void initialise(const int nthreads = 0);
  void exit();
  int get_number_of_processors();
  int get_processor_id();

}

#endif // MPI_WRAPPER
