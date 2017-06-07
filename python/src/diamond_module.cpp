
//#include <boost/python.hpp>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
//#include "include/numpy_boost_python.hpp"

#include "diamond_wrapper.hpp"
using namespace boost::python;

void export_reconstruction();
void export_filters();
BOOST_PYTHON_MODULE(parallelbeam)
{
	namespace bp = boost::python;
	np::initialize();
	//To specify that this module is a package
	bp::object package = bp::scope();
	package.attr("__path__") = "parallelbeam";

	export_reconstruction();
	export_filters();
}

void export_reconstruction()
{
	namespace bp = boost::python;
	bp::object reconstructionModule(bp::handle<>(bp::borrowed(PyImport_AddModule("parallelbeam.alg"))));
	// make "from diamond import filters" work
	bp::scope().attr("alg") = reconstructionModule;	
	// set the current scope to the new sub-module
    bp::scope reconstruction_scope = reconstructionModule;

  def("cgls", reconstruct_cgls);
  def("cgls_conv", reconstruct_cgls2);
  def("sirt", reconstruct_sirt);
  def("mlem", reconstruct_mlem);
  def("cgls_tikhonov", reconstruct_cgls_tikhonov);
  def("cgls_TVreg", reconstruct_cgls_tvreg);
  //def("tvreg", reconstruct_tvreg);
}

void export_filters()
{
	namespace bp = boost::python;
	bp::object filtersModule(bp::handle<>(bp::borrowed(PyImport_AddModule("parallelbeam.filters"))));
	// make "from diamond import filters" work
	bp::scope().attr("filters") = filtersModule;	
	// set the current scope to the new sub-module
    bp::scope filters_scope = filtersModule;

	def("aml_ring_artefacts", ring_artefacts_aml);
}