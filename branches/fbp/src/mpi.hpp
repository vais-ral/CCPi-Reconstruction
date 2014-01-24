
#ifndef MPI_WRAPPER
#define MPI_WRAPPER

namespace machine {

  void initialise();
  void exit();
  int get_number_of_processors();
  int get_processor_id();

}

#endif // MPI_WRAPPER
