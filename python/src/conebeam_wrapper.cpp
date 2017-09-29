#include <iostream>
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


np::ndarray conebeam_reconstruct_iter(np::ndarray ndarray_pixels,
				       np::ndarray ndarray_angles,
					   np::ndarray ndarray_h_offsets,
					   np::ndarray ndarray_v_offsets,
					   const int pixels_per_voxel, const double source_x,
					   const double detector_x, const double pixel_h_size,
					   const double pixel_v_size, const double mask_radius,
					   const bool beam_harden,
					   np::ndarray ndarray_full_vox_origin,
					   np::ndarray ndarray_voxel_size,
				       int niterations, int nthreads,
				       CCPi::algorithms alg,					   
				       const real regularize = 0.0,
					   np::ndarray norm=np::zeros(bp::make_tuple(0), np::dtype::get_builtin<float>()),					   
					   bool is_pixels_in_log=true)
{
  bool beam_hardening = false; //overwriting as beam hardening is not implemented
  // vertical size to break data up into for processing
  const int blocking_factor = 0;
  // number of GPUs etc if using accelerated code
  //const int num_devices = 1;
  numpy_3d pixels(reinterpret_cast<float*>(ndarray_pixels.get_data()), boost::extents[ndarray_pixels.shape(0)][ndarray_pixels.shape(1)][ndarray_pixels.shape(2)]);  
  numpy_1d angles(reinterpret_cast<float*>(ndarray_angles.get_data()), boost::extents[ndarray_angles.shape(0)]);  
  numpy_1d h_offsets(reinterpret_cast<float*>(ndarray_h_offsets.get_data()), boost::extents[ndarray_h_offsets.shape(0)]);  
  numpy_1d v_offsets(reinterpret_cast<float*>(ndarray_v_offsets.get_data()), boost::extents[ndarray_v_offsets.shape(0)]);    
  numpy_1d full_vox_origin(reinterpret_cast<float*>(ndarray_full_vox_origin.get_data()), boost::extents[ndarray_full_vox_origin.shape(0)]);    
  numpy_1d voxel_size(reinterpret_cast<float*>(ndarray_voxel_size.get_data()), boost::extents[ndarray_voxel_size.shape(0)]);      				   
					   
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
	  norm[i] = 0.0;
      } else {
	real_1d data(niterations);
	for (int i = 0; i < niterations; i++)
	  norm[i] = data[i];
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
  np::dtype dt = np::dtype::get_builtin<float>();
  bp::tuple shape = bp::make_tuple(dims[0], dims[1], dims[2]);
  np::ndarray varray = np::zeros(shape,dt);
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


np::ndarray conebeam_reconstruct_cgls(np::ndarray pixels,
				       np::ndarray angles,
					   np::ndarray h_offsets,
					   np::ndarray v_offsets,
					   const int pixels_per_voxel, const double source_x,
					   const double detector_x, const double pixel_h_size,
					   const double pixel_v_size, const double mask_radius,
					   const bool beam_harden, np::ndarray full_vox_origin,
					   np::ndarray voxel_size,
				       int niterations, int nthreads,
					   bool is_pixels_in_log=true)
{
  return conebeam_reconstruct_iter(pixels, angles, h_offsets, v_offsets, pixels_per_voxel, source_x, detector_x, pixel_h_size, pixel_v_size,mask_radius,
			  beam_harden, full_vox_origin, voxel_size, niterations, nthreads, CCPi::alg_CGLS,  0.0, np::zeros(bp::make_tuple(1), np::dtype::get_builtin<float>()), is_pixels_in_log);
}

np::ndarray conebeam_reconstruct_sirt(np::ndarray pixels,
				       np::ndarray angles,
					   np::ndarray h_offsets,
					   np::ndarray v_offsets,
					   const int pixels_per_voxel, const double source_x,
					   const double detector_x, const double pixel_h_size,
					   const double pixel_v_size, const double mask_radius,
					   const bool beam_harden, np::ndarray full_vox_origin,
					   np::ndarray voxel_size,
				       int niterations, int nthreads,
					   bool is_pixels_in_log=true)
{
  return conebeam_reconstruct_iter(pixels, angles, h_offsets, v_offsets, pixels_per_voxel, source_x, detector_x, pixel_h_size, pixel_v_size,mask_radius,
			  beam_harden, full_vox_origin, voxel_size, niterations, nthreads, CCPi::alg_SIRT,  0.0, np::zeros(bp::make_tuple(1), np::dtype::get_builtin<float>()), is_pixels_in_log);
}

np::ndarray conebeam_reconstruct_mlem(np::ndarray pixels,
				       np::ndarray angles,
					   np::ndarray h_offsets,
					   np::ndarray v_offsets,
					   const int pixels_per_voxel, const double source_x,
					   const double detector_x, const double pixel_h_size,
					   const double pixel_v_size, const double mask_radius,
					   const bool beam_harden, np::ndarray full_vox_origin,
					   np::ndarray voxel_size,
				       int niterations, int nthreads,
					   bool is_pixels_in_log=true)
{
  return conebeam_reconstruct_iter(pixels, angles, h_offsets, v_offsets, pixels_per_voxel, source_x, detector_x, pixel_h_size, pixel_v_size,mask_radius,
			  beam_harden, full_vox_origin, voxel_size, niterations, nthreads, CCPi::alg_MLEM, 0.0, np::zeros(bp::make_tuple(1), np::dtype::get_builtin<float>()),  is_pixels_in_log);	
}

np::ndarray
conebeam_reconstruct_cgls_tikhonov(np::ndarray pixels,
				       np::ndarray angles,
					   np::ndarray h_offsets,
					   np::ndarray v_offsets,
					   const int pixels_per_voxel, const double source_x,
					   const double detector_x, const double pixel_h_size,
					   const double pixel_v_size, const double mask_radius,
					   const bool beam_harden, np::ndarray full_vox_origin,
					   np::ndarray voxel_size,
				       int niterations, int nthreads,
					   double regularize,
					   np::ndarray norm_r, 
					   bool is_pixels_in_log=true)
{
  return conebeam_reconstruct_iter(pixels, angles, h_offsets, v_offsets, pixels_per_voxel, source_x, detector_x, pixel_h_size, pixel_v_size,mask_radius,
			  beam_harden, full_vox_origin, voxel_size, niterations, nthreads, CCPi::alg_CGLS_Tikhonov,  regularize, norm_r,  is_pixels_in_log);	
}

np::ndarray
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
					   bool is_pixels_in_log=true)
{
  return conebeam_reconstruct_iter(pixels, angles, h_offsets, v_offsets, pixels_per_voxel, source_x, detector_x, pixel_h_size, pixel_v_size,mask_radius,
			  beam_harden, full_vox_origin, voxel_size, niterations, nthreads, CCPi::alg_CGLS_TVreg,  regularize,  norm_r, is_pixels_in_log);		
}

np::ndarray conebeam_reconstruct_cgls2(np::ndarray pixels,
				       np::ndarray angles,
					   np::ndarray h_offsets,
					   np::ndarray v_offsets,
					   const int pixels_per_voxel, const double source_x,
					   const double detector_x, const double pixel_h_size,
					   const double pixel_v_size, const double mask_radius,
					   const bool beam_harden, np::ndarray full_vox_origin,
					   np::ndarray voxel_size,
				       int niterations, int nthreads,
					   np::ndarray norm_r, 
					   bool is_pixels_in_log=true)
{
  return conebeam_reconstruct_iter(pixels, angles, h_offsets, v_offsets, pixels_per_voxel, source_x, detector_x, pixel_h_size, pixel_v_size,mask_radius,
			  beam_harden, full_vox_origin, voxel_size, niterations, nthreads, CCPi::alg_CGLS,  0.0, norm_r,
			   is_pixels_in_log);	
}

void reconstruct_tvreg()
{
}


boost::python::tuple conebeam_create_phantom()
{
	CCPi::Nikon_XTek *instrument = new CCPi::Nikon_XTek();
	instrument->initialise_phantom();
	int nx = instrument->get_num_angles();
	int ny = instrument->get_num_h_pixels();
	int nz = instrument->get_num_v_pixels();
	double source_x = -1*instrument->get_source_x();
	double detector_x = instrument->get_detector_x()+source_x;
	double mask_radius = instrument->get_mask_radius();
	double pixel_h_size = 0.390625;
	double pixel_v_size = 0.390625;
	std::cout<<" dimension "<<nx<<" "<<ny<<" "<<std::endl;

	np::dtype dt = np::dtype::get_builtin<float>();
	bp::tuple shape = bp::make_tuple(nx, ny, nz);
	np::ndarray pixels = np::zeros(shape,dt);
	
	pixel_3d data = instrument->get_pixel_data();
	for(int i=0;i<nx;i++)
		for(int j=0;j<ny;j++)
			for(int k=0;k<nz;k++)
				pixels[i][j][k] = data[i][j][k];

	bp::tuple ashape = bp::make_tuple(instrument->get_num_angles());
	np::ndarray angles = np::zeros(ashape,dt);	
	real_1d phi = instrument->get_phi();
	for(int i=0;i<instrument->get_num_angles();i++)
		angles[i] = phi[i]*180.0/M_PI;
	return boost::python::make_tuple(pixels, angles, source_x, detector_x, pixel_h_size, pixel_v_size, mask_radius);
}