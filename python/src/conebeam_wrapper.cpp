#include "conebeam_wrapper.hpp"
#include "algorithms.hpp"
#include "instruments.hpp"
#include "cgls.hpp"
#include "mlem.hpp"
#include "sirt.hpp"
#include "tv_reg.hpp"
#include "results.hpp"
#include "voxels.hpp"
#include "mpi.hpp"


//voxel_data *reconstruct(CCPi::instrument *device,
//			CCPi::reconstruction_alg *algorithm,
//			const numpy_3d &pixels, const numpy_1d &angles,
//			const numpy_1d &h_offsets, const numpy_1d &v_offsets,
//			const int pixels_per_voxel, const real source_x,
//			const real detector_x, const real pixel_h_size,
//			const real pixel_v_size, const real mask_radius,
//			const bool beam_harden, real full_vox_origin[3],
//			real voxel_size[3], const bool has_offsets,bool is_pixels_in_log)

numpy_boost<float, 3> conebeam_reconstruct_iter(const numpy_boost<float, 3> &pixels,
				       const numpy_boost<float, 1> &angles,
					   const numpy_boost<float, 1> &h_offsets,
					   const numpy_boost<float, 1> &v_offsets,
					   const int pixels_per_voxel, const double source_x,
					   const double detector_x, const double pixel_h_size,
					   const double pixel_v_size, const double mask_radius,
					   const bool beam_harden,
					   const numpy_boost<float, 1> &full_vox_origin,
					   const numpy_boost<float, 1> &voxel_size,
				       int niterations, int nthreads,
				       CCPi::algorithms alg,
				       const real regularize = 0.0,
					   numpy_boost<float, 1> *norm = 0, bool is_pixels_in_log=true)
{
  bool beam_hardening = false; //overwriting as beam hardening is not implemented
  // vertical size to break data up into for processing
  const int blocking_factor = 0;
  // number of GPUs etc if using accelerated code
  //const int num_devices = 1;
  voxel_data *voxels = 0;
  CCPi::reconstruction_alg *algorithm = 0;
  switch (alg) {
  case CCPi::alg_CGLS:
    //if (blocking_factor > 0 and instrument->supports_blocks())
    //  recon_algorithm = new CCPi::cgls_2d(niterations, pixels_per_voxel);
    algorithm = new CCPi::cgls_3d(niterations);
    break;
  case CCPi::alg_SIRT:
    algorithm = new CCPi::sirt(niterations);
    break;
  case CCPi::alg_MLEM:
    algorithm = new CCPi::mlem(niterations);
    break;
  case CCPi::alg_CGLS_Tikhonov:
    algorithm = new CCPi::cgls_tikhonov(niterations, regularize);
    break;
  case CCPi::alg_CGLS_TVreg:
    algorithm = new CCPi::cgls_tv_reg(niterations, regularize);
    break;
  default:
    break;
  }
  if (algorithm != 0) {
	CCPi::instrument *instrument = new CCPi::Nikon_XTek();
    machine::initialise(nthreads);
	real vox_origin[3];
	real vox_size[3];
	vox_origin[0]=full_vox_origin[0];vox_origin[1]=full_vox_origin[1];vox_origin[2]=full_vox_origin[2]; 
	vox_size[0]=voxel_size[0];vox_size[1]=voxel_size[1];vox_size[2]=voxel_size[2];
	
    // instrument setup from pixels/angles will probably copy
	voxels = reconstruct(instrument, algorithm, pixels,
					 angles, h_offsets, v_offsets,
					 pixels_per_voxel, source_x, detector_x,
					 pixel_h_size, pixel_v_size, mask_radius,
					 beam_hardening, vox_origin, vox_size,
					 false, is_pixels_in_log);	
    //voxels = reconstruct(instrument, algorithm, pixels, angles, rotation_centre,
	//		 resolution, blocking_factor, beam_harden,is_pixels_in_log);
    machine::exit();
    delete instrument;
    if (norm != 0) {
      if (voxels == 0) {
	for (int i = 0; i < niterations; i++)
	  (*norm)[i] = 0.0;
      } else {
	real_1d data(niterations);
	for (int i = 0; i < niterations; i++)
	  (*norm)[i] = data[i];
      }
    }
    delete algorithm;
  }
  int dims[3];
  if (voxels == 0) {
    dims[0] = 1;
    dims[1] = 1;
    dims[2] = 1;
  } else {
    // Todo - remove buffered region
    dims[0] = voxels->shape()[0];
    dims[1] = voxels->shape()[1];
    dims[2] = voxels->shape()[2];
  }
  numpy_boost<float, 3> varray(dims);
  if (voxels == 0)
    varray[0][0][0] = 0.0;
  else {
    // Todo - vector? and remove buffered region
    for (int i = 0; i < dims[0]; i++)
      for (int j = 0; j < dims[1]; j++)
	for (int k = 0; k < dims[2]; k++)
	  varray[i][j][k] = (*voxels)[i][j][k];
    delete voxels;
  }
  return varray;
}


numpy_boost<float, 3> conebeam_reconstruct_cgls(const numpy_boost<float, 3> &pixels,
				       const numpy_boost<float, 1> &angles,
					   const numpy_boost<float, 1> &h_offsets,
					   const numpy_boost<float, 1> &v_offsets,
					   const int pixels_per_voxel, const double source_x,
					   const double detector_x, const double pixel_h_size,
					   const double pixel_v_size, const double mask_radius,
					   const bool beam_harden, const numpy_boost<float, 1> &full_vox_origin,
					   const numpy_boost<float, 1> &voxel_size,
				       int niterations, int nthreads,
					   bool is_pixels_in_log=true)
{
  return conebeam_reconstruct_iter(pixels, angles, h_offsets, v_offsets, pixels_per_voxel, source_x, detector_x, pixel_h_size, pixel_v_size,mask_radius,
			  beam_harden, full_vox_origin, voxel_size, niterations, nthreads, CCPi::alg_CGLS, 0.0, 0, is_pixels_in_log);
}

numpy_boost<float, 3> conebeam_reconstruct_sirt(const numpy_boost<float, 3> &pixels,
				       const numpy_boost<float, 1> &angles,
					   const numpy_boost<float, 1> &h_offsets,
					   const numpy_boost<float, 1> &v_offsets,
					   const int pixels_per_voxel, const double source_x,
					   const double detector_x, const double pixel_h_size,
					   const double pixel_v_size, const double mask_radius,
					   const bool beam_harden, const numpy_boost<float, 1> &full_vox_origin,
					   const numpy_boost<float, 1> &voxel_size,
				       int niterations, int nthreads,
					   bool is_pixels_in_log=true)
{
  return conebeam_reconstruct_iter(pixels, angles, h_offsets, v_offsets, pixels_per_voxel, source_x, detector_x, pixel_h_size, pixel_v_size,mask_radius,
			  beam_harden, full_vox_origin, voxel_size, niterations, nthreads, CCPi::alg_SIRT, 0.0, 0, is_pixels_in_log);
}

numpy_boost<float, 3> conebeam_reconstruct_mlem(const numpy_boost<float, 3> &pixels,
				       const numpy_boost<float, 1> &angles,
					   const numpy_boost<float, 1> &h_offsets,
					   const numpy_boost<float, 1> &v_offsets,
					   const int pixels_per_voxel, const double source_x,
					   const double detector_x, const double pixel_h_size,
					   const double pixel_v_size, const double mask_radius,
					   const bool beam_harden, const numpy_boost<float, 1> &full_vox_origin,
					   const numpy_boost<float, 1> &voxel_size,
				       int niterations, int nthreads,
					   bool is_pixels_in_log=true)
{
  return conebeam_reconstruct_iter(pixels, angles, h_offsets, v_offsets, pixels_per_voxel, source_x, detector_x, pixel_h_size, pixel_v_size,mask_radius,
			  beam_harden, full_vox_origin, voxel_size, niterations, nthreads, CCPi::alg_MLEM, 0.0, 0, is_pixels_in_log);	
}

numpy_boost<float, 3>
conebeam_reconstruct_cgls_tikhonov(const numpy_boost<float, 3> &pixels,
				       const numpy_boost<float, 1> &angles,
					   const numpy_boost<float, 1> &h_offsets,
					   const numpy_boost<float, 1> &v_offsets,
					   const int pixels_per_voxel, const double source_x,
					   const double detector_x, const double pixel_h_size,
					   const double pixel_v_size, const double mask_radius,
					   const bool beam_harden, const numpy_boost<float, 1> &full_vox_origin,
					   const numpy_boost<float, 1> &voxel_size,
				       int niterations, int nthreads,
					   double regularize,
					   numpy_boost<float, 1> norm_r, 
					   bool is_pixels_in_log=true)
{
  return conebeam_reconstruct_iter(pixels, angles, h_offsets, v_offsets, pixels_per_voxel, source_x, detector_x, pixel_h_size, pixel_v_size,mask_radius,
			  beam_harden, full_vox_origin, voxel_size, niterations, nthreads, CCPi::alg_CGLS_Tikhonov, regularize, &norm_r, is_pixels_in_log);	
}

numpy_boost<float, 3>
conebeam_reconstruct_cgls_tvreg(const numpy_boost<float, 3> &pixels,
				       const numpy_boost<float, 1> &angles,
					   const numpy_boost<float, 1> &h_offsets,
					   const numpy_boost<float, 1> &v_offsets,
					   const int pixels_per_voxel, 
					   const double source_x,
					   const double detector_x, 
					   const double pixel_h_size,
					   const double pixel_v_size, 
					   const double mask_radius,
					   const bool beam_harden, 
					   const numpy_boost<float, 1> &full_vox_origin,
					   const numpy_boost<float, 1> &voxel_size,
				       int niterations, 
					   int nthreads,
					   double regularize,
					   numpy_boost<float, 1> norm_r, 
					   bool is_pixels_in_log=true)
{
  return conebeam_reconstruct_iter(pixels, angles, h_offsets, v_offsets, pixels_per_voxel, source_x, detector_x, pixel_h_size, pixel_v_size,mask_radius,
			  beam_harden, full_vox_origin, voxel_size, niterations, nthreads, CCPi::alg_CGLS_TVreg, regularize, &norm_r, is_pixels_in_log);		
}

numpy_boost<float, 3> conebeam_reconstruct_cgls2(const numpy_boost<float, 3> &pixels,
				       const numpy_boost<float, 1> &angles,
					   const numpy_boost<float, 1> &h_offsets,
					   const numpy_boost<float, 1> &v_offsets,
					   const int pixels_per_voxel, const double source_x,
					   const double detector_x, const double pixel_h_size,
					   const double pixel_v_size, const double mask_radius,
					   const bool beam_harden, const numpy_boost<float, 1> &full_vox_origin,
					   const numpy_boost<float, 1> &voxel_size,
				       int niterations, int nthreads,
					   numpy_boost<float, 1> norm_r, 
					   bool is_pixels_in_log=true)
{
  return conebeam_reconstruct_iter(pixels, angles, h_offsets, v_offsets, pixels_per_voxel, source_x, detector_x, pixel_h_size, pixel_v_size,mask_radius,
			  beam_harden, full_vox_origin, voxel_size, niterations, nthreads, CCPi::alg_CGLS, 0.0,
			  &norm_r, is_pixels_in_log);	
}

void reconstruct_tvreg()
{
}
