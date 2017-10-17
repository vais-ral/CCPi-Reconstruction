
//#include <boost/python.hpp>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
//This will enable more than 15 arguments in function
#define BOOST_PYTHON_MAX_ARITY 19  

#include "conebeam_wrapper.hpp"
using namespace boost::python;

void export_conebeam_reconstruction();

BOOST_PYTHON_MODULE(conebeam)
{
	namespace bp = boost::python;
	np::initialize();	
	//To specify that this module is a package
	bp::object package = bp::scope();
	package.attr("__path__") = "conebeam";
			
	export_conebeam_reconstruction();
}

void export_conebeam_reconstruction()
{
	namespace bp = boost::python;
	bp::object reconstructionModule(bp::handle<>(bp::borrowed(PyImport_AddModule("conebeam.alg"))));
	// make "from diamond import filters" work
	bp::scope().attr("alg") = reconstructionModule;	
	// set the current scope to the new sub-module
    bp::scope reconstruction_scope = reconstructionModule;

  def("cgls", conebeam_reconstruct_cgls);
  def("cgls_conv", conebeam_reconstruct_cgls2);
  def("sirt", conebeam_reconstruct_sirt);
  def("mlem", conebeam_reconstruct_mlem);
  def("cgls_tikhonov", conebeam_reconstruct_cgls_tikhonov);
  def("cgls_TVreg", conebeam_reconstruct_cgls_tvreg);
  def("create_phantom", conebeam_create_phantom);
  def("load_xtek", conebeam_load_xtek);
}
