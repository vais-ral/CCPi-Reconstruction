#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include "base_types.hpp"

namespace bp = boost::python;
namespace np = boost::python::numpy;
extern np::ndarray
ring_artefacts_aml(np::ndarray pixels,
		   const float param_n, const float param_r,
		   const int num_series);
extern np::ndarray
reconstruct_cgls(np::ndarray pixels,
		 np::ndarray angles,
		 double rotation_centre, int resolution,
		 int niterations, int nthreads,bool is_pixels_in_log);

extern np::ndarray
reconstruct_cgls_step(np::ndarray pixels,
	np::ndarray angles,
	double rotation_centre, int resolution,
	int niterations, int nthreads, bool is_pixels_in_log,
	np::ndarray last_iteration_volume);

extern np::ndarray
reconstruct_sirt(np::ndarray pixels,
		 np::ndarray angles,
		 double rotation_centre, int resolution,
		 int niterations, int nthreads,bool is_pixels_in_log);
extern np::ndarray
reconstruct_sirt_step(np::ndarray pixels,
	np::ndarray angles,
	double rotation_centre, int resolution,
	int niterations, int nthreads, bool is_pixels_in_log,
	np::ndarray last_iteration_volume);


extern np::ndarray
reconstruct_mlem(np::ndarray pixels,
		 np::ndarray angles,
		 double rotation_centre, int resolution,
		 int niterations, int nthreads,bool is_pixels_in_log);
extern np::ndarray
reconstruct_mlem_step(np::ndarray pixels,
	np::ndarray angles,
	double rotation_centre, int resolution,
	int niterations, int nthreads, bool is_pixels_in_log,
	np::ndarray last_iteration_volume);

extern np::ndarray
reconstruct_cgls2(np::ndarray pixels,
		  np::ndarray angles,
		  double rotation_centre, int resolution,
		  int niterations, int nthreads, np::ndarray norm_r,bool is_pixels_in_log);
extern np::ndarray
reconstruct_cgls2_step(np::ndarray pixels,
	np::ndarray angles,
	double rotation_centre, int resolution,
	int niterations, int nthreads, np::ndarray norm_r, bool is_pixels_in_log,
	np::ndarray last_iteration_volume);

extern np::ndarray
reconstruct_cgls_tikhonov(np::ndarray pixels,
			  np::ndarray angles,
			  double rotation_centre, int resolution,
			  int niterations, int nthreads, double regularize,
			  np::ndarray norm_r,bool is_pixels_in_log);
extern np::ndarray
reconstruct_cgls_tikhonov_step(np::ndarray pixels,
	np::ndarray angles,
	double rotation_centre, int resolution,
	int niterations, int nthreads, double regularize,
	np::ndarray norm_r, bool is_pixels_in_log,
	np::ndarray last_iteration_volume);

extern np::ndarray
reconstruct_cgls_tvreg(np::ndarray pixels,
		       np::ndarray angles,
		       double rotation_centre, int resolution,
		       int niterations, int nthreads, double regularize,
		       np::ndarray norm_r,bool is_pixels_in_log);
extern np::ndarray
reconstruct_cgls_tvreg_step(np::ndarray pixels,
	np::ndarray angles,
	double rotation_centre, int resolution,
	int niterations, int nthreads, double regularize,
	np::ndarray norm_r, bool is_pixels_in_log,
	np::ndarray last_iteration_volume);

extern void reconstruct_tvreg();
