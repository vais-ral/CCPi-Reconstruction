
//#include <boost/python.hpp>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "include/numpy_boost_python.hpp"

#include "base_types.hpp"
#include "mpi.hpp"
#include "utils.hpp"
#include "instruments.hpp"
#include "algorithms.hpp"
#include "results.hpp"
#include "voxels.hpp"

using namespace boost::python;

#include "reconstruct.hpp"
#include "reconstruct.cpp"

BOOST_PYTHON_MODULE(ccpi)
{
  import_array();
  numpy_boost_python_register_type<float, 1>();
  //numpy_boost_python_register_type<float, 2>();
  numpy_boost_python_register_type<float, 3>();
  //def("test", test);
  def("cgls", reconstruct_cgls);
  def("tvreg", reconstruct_tvreg);
}
