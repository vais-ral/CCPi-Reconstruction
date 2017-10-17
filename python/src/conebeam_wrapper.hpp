#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include "base_types.hpp"

namespace bp = boost::python;
namespace np = boost::python::numpy;

extern np::ndarray conebeam_reconstruct_cgls(np::ndarray pixels,
				       np::ndarray angles,
					   np::ndarray h_offsets,
					   np::ndarray v_offsets,
					   const int pixels_per_voxel,
					   const double source_x,
					   const double detector_x,
					   const double pixel_h_size,
					   const double pixel_v_size,
					   const double mask_radius,
					   const bool beam_harden,
					   np::ndarray full_vox_origin,
					   np::ndarray voxel_size,
				       int niterations, 
					   int nthreads,
					   bool is_pixels_in_log);
					   
extern np::ndarray conebeam_reconstruct_sirt(np::ndarray pixels,
				       np::ndarray angles,
					   np::ndarray h_offsets,
					   np::ndarray v_offsets,
					   const int pixels_per_voxel,
					   const double source_x,
					   const double detector_x, 
					   const double pixel_h_size,
					   const double pixel_v_size, 
					   const double mask_radius,
					   const bool beam_harden, 
					   np::ndarray full_vox_origin,
					   np::ndarray voxel_size,
				       int niterations, int nthreads,
					   bool is_pixels_in_log);
					   
extern np::ndarray conebeam_reconstruct_mlem(np::ndarray pixels,
				       np::ndarray angles,
					   np::ndarray h_offsets,
					   np::ndarray v_offsets,
					   const int pixels_per_voxel, const double source_x,
					   const double detector_x, const double pixel_h_size,
					   const double pixel_v_size, const double mask_radius,
					   const bool beam_harden, 
					   np::ndarray full_vox_origin,
					   np::ndarray voxel_size,
				       int niterations, int nthreads,
					   bool is_pixels_in_log);
					   
extern np::ndarray
conebeam_reconstruct_cgls_tikhonov(np::ndarray pixels,
				       np::ndarray angles,
					   np::ndarray h_offsets,
					   np::ndarray v_offsets,
					   const int pixels_per_voxel, 
					   const double source_x,
					   const double detector_x, 
					   const double pixel_h_size,
					   const double pixel_v_size, 
					   const double mask_radius,
					   const bool beam_harden, 
					   np::ndarray full_vox_origin,
					   np::ndarray voxel_size,
				       int niterations, 
					   int nthreads,
					   double regularize,
					   np::ndarray norm_r, 
					   bool is_pixels_in_log);
					   
extern np::ndarray
conebeam_reconstruct_cgls_tvreg(np::ndarray pixels,
				       np::ndarray angles,
					   np::ndarray h_offsets,
					   np::ndarray v_offsets,
					   const int pixels_per_voxel, 
					   const double source_x,
					   const double detector_x, 
					   const double pixel_h_size,
					   const double pixel_v_size, 
					   const double mask_radius,
					   const bool beam_harden, 
					   np::ndarray full_vox_origin,
					   np::ndarray voxel_size,
				       int niterations, 
					   int nthreads,
					   double regularize,
					   np::ndarray norm_r, 
					   bool is_pixels_in_log);
extern np::ndarray conebeam_reconstruct_cgls2(np::ndarray pixels,
				       np::ndarray angles,
					   np::ndarray h_offsets,
					   np::ndarray v_offsets,
					   const int pixels_per_voxel, 
					   const double source_to_sample,
					   const double source_to_detector, 
					   const double pixel_h_size,
					   const double pixel_v_size, 
					   const double mask_radius,
					   const bool beam_harden, 
					   np::ndarray full_vox_origin,
					   np::ndarray voxel_size,
				       int niterations, int nthreads,
					   np::ndarray norm_r, 
					   bool is_pixels_in_log);
					   
extern boost::python::tuple conebeam_create_phantom();
extern boost::python::tuple conebeam_load_xtek(const std::string& filename);
