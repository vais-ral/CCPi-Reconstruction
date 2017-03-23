
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
#include "cgls.hpp"
#include "sirt.hpp"
#include "mlem.hpp"

using namespace boost::python;

#include "reconstruct.hpp"
#include "reconstruct.cpp"

BOOST_PYTHON_MODULE(reconstruction)
{
  import_array();
  numpy_boost_python_register_type<float, 1>();
  //numpy_boost_python_register_type<float, 2>();
  numpy_boost_python_register_type<float, 3>();
  numpy_boost_python_register_type<double, 3>();
  def("cgls", reconstruct_cgls);
  def("cgls_conv", reconstruct_cgls2);
  def("sirt", reconstruct_sirt);
  def("mlem", reconstruct_mlem);
  def("cgls_tikhonov", reconstruct_cgls_tikhonov);
  def("cgls_TVreg", reconstruct_cgls_tvreg);
  //def("tvreg", reconstruct_tvreg);
}

