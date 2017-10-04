#include "diamond_wrapper.hpp"
#include "algorithms.hpp"
#include "instruments.hpp"
#include "cgls.hpp"
#include "mlem.hpp"
#include "sirt.hpp"
#include "tv_reg.hpp"
#include "results.hpp"
#include "voxels.hpp"
#include "mpi.hpp"
#include "base_types.hpp"
#include "utils.hpp"

// assuming normalised data from SAVU
// Inputs : numpy_boost<double,3> pixels
// Returns : double,3
np::ndarray ring_artefacts_aml(np::ndarray pixels,
					  const float param_n,
					  const float param_r,
					  const int num_series)
{
  // Todo - worry about float/double precision loss
  sl_int nangles = (sl_int)pixels.shape(0);
  sl_int nh = (sl_int)pixels.shape(2);
  sl_int nv = (sl_int)pixels.shape(1);
  // issue with change of precision?
  pixel_3d p(boost::extents[nangles][nh][nv]);
  for (int i = 0; i < nangles; i++) {
	for (sl_int j = 0; j < nh; j++) {
		for (sl_int k = 0; k < nv; k++) {
			p[i][j][k] = - std::log(bp::extract<float>(pixels[i][k][j]));
		}
	}
  }	  
  CCPi::remove_aml_ring_artefacts(p, nangles, nh, nv, param_n,
				  param_r, num_series);
  int dims[3];
  dims[0] = nangles;
  dims[1] = nv;
  dims[2] = nh;
  np::dtype dt = np::dtype::get_builtin<float>();
  bp::tuple shape = bp::make_tuple(dims[0], dims[1], dims[2]); 
  np::ndarray varray = np::zeros(shape,dt);
  // Todo - vector? and remove buffered region
  for (int i = 0; i < dims[0]; i++)
    for (int j = 0; j < dims[1]; j++)
      for (int k = 0; k < dims[2]; k++)
	varray[i][j][k] = std::exp(- p[i][k][j]);
  return varray;
}

np::ndarray reconstruct_iter(np::ndarray ndarray_pixels,
				       np::ndarray ndarray_angles,
				       double rotation_centre, int resolution,
				       int niterations, int nthreads,
				       CCPi::algorithms alg,
				       const real regularize = 0.0,
				       np::ndarray norm=np::zeros(bp::make_tuple(0), np::dtype::get_builtin<float>()),bool is_pixels_in_log=true)
{
  bool beam_harden = false;
  // vertical size to break data up into for processing
  const int blocking_factor = 0;
  // number of GPUs etc if using accelerated code
  //const int num_devices = 1;
  numpy_3d pixels(reinterpret_cast<float*>(ndarray_pixels.get_data()), boost::extents[ndarray_pixels.shape(0)][ndarray_pixels.shape(1)][ndarray_pixels.shape(2)]);  
  numpy_1d angles(reinterpret_cast<float*>(ndarray_angles.get_data()), boost::extents[ndarray_angles.shape(0)]);
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
    CCPi::instrument *instrument = new CCPi::Diamond();
    machine::initialise(nthreads);
    // instrument setup from pixels/angles will probably copy
    voxels = reconstruct(instrument, algorithm, pixels, angles, rotation_centre,
			 resolution, blocking_factor, beam_harden,is_pixels_in_log);
    machine::exit();
    delete instrument;
    if (norm[0] != 0) {
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
/*********WITH PREVIOUS STEP************/
np::ndarray reconstruct_iter(np::ndarray ndarray_pixels,
	np::ndarray ndarray_angles,
	double rotation_centre, int resolution,
	int niterations, int nthreads,
	CCPi::algorithms alg,
	np::ndarray last_iteration_volume ,
	const real regularize = 0.0,
	np::ndarray norm = np::zeros(bp::make_tuple(0), np::dtype::get_builtin<float>()), bool is_pixels_in_log = true
	)
{
	bool beam_harden = false;
	// vertical size to break data up into for processing
	const int blocking_factor = 0;
	// number of GPUs etc if using accelerated code
	//const int num_devices = 1;
	numpy_3d pixels(reinterpret_cast<float*>(ndarray_pixels.get_data()), boost::extents[ndarray_pixels.shape(0)][ndarray_pixels.shape(1)][ndarray_pixels.shape(2)]);
	numpy_1d angles(reinterpret_cast<float*>(ndarray_angles.get_data()), boost::extents[ndarray_angles.shape(0)]);
	voxel_data *voxels = 0;
	/******************************************************************************************/
	/*
	Copy the last_iteration_volume numpy array into the
	voxels boost multi array
	*/
	int nx_voxels = last_iteration_volume.shape(0);
	int ny_voxels = last_iteration_volume.shape(1);
	int nz_voxels = last_iteration_volume.shape(2);
	
	printf("last_iteration_volume size [%d,%d,%d]\n",nx_voxels, ny_voxels, nz_voxels );
	voxel_type * A = reinterpret_cast<voxel_type *>(last_iteration_volume.get_data());

	voxel_data * input_voxels = new voxel_data(boost::extents[nx_voxels][ny_voxels][nz_voxels]);

	for (int i = 0; i < nx_voxels; i++) {
		for (int j = 0; j < ny_voxels; j++) {
			for (int k = 0; k < nz_voxels; k++) {
				int index = k + nz_voxels * j + nz_voxels * ny_voxels * i;
				voxel_type val = (*(A + index));
				(*input_voxels)[i][j][k] = val;
				//std::memcpy((*input_voxels)[i][j][k], &val, sizeof(voxel_type));
			}
		}
	}


	/******************************************************************************************/
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
		CCPi::instrument *instrument = new CCPi::Diamond();
		
		machine::initialise(nthreads);
		// instrument setup from pixels/angles will probably copy
		voxels = reconstruct(instrument, algorithm, pixels, angles, rotation_centre,
			resolution, blocking_factor, beam_harden, is_pixels_in_log, input_voxels);
		machine::exit();
		delete instrument;
		if (norm[0] != 0) {
			if (voxels == 0) {
				for (int i = 0; i < niterations; i++)
					norm[i] = 0.0;
			}
			else {
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
	}
	else {
		// Todo - remove buffered region
		dims[0] = voxels->shape()[0];
		dims[1] = voxels->shape()[1];
		dims[2] = voxels->shape()[2];
	}

	np::dtype dt = np::dtype::get_builtin<float>();
	bp::tuple shape = bp::make_tuple(dims[0], dims[1], dims[2]);
	np::ndarray varray = np::zeros(shape, dt);
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
	//QUESTION: this line makes python crash. Possibly because it destroys memory from input numpy array?
	//delete input_voxels;
	return varray;
}

np::ndarray reconstruct_cgls(np::ndarray pixels,
				       np::ndarray angles,
				       double rotation_centre, int resolution,
				       int niterations, int nthreads,bool is_pixels_in_log)
{
  return reconstruct_iter(pixels, angles, rotation_centre, resolution,
			  niterations, nthreads, CCPi::alg_CGLS, 0.0, np::zeros(bp::make_tuple(1), np::dtype::get_builtin<float>()), is_pixels_in_log);
}
/*If you provide a previous step*/
np::ndarray reconstruct_cgls_step(np::ndarray pixels,
	np::ndarray angles,
	double rotation_centre, int resolution,
	int niterations, int nthreads, bool is_pixels_in_log,
	np::ndarray last_iteration_volume)
{
	return reconstruct_iter(pixels, angles, rotation_centre, resolution,
		niterations, nthreads, CCPi::alg_CGLS, last_iteration_volume, 
		0.0, np::zeros(bp::make_tuple(1), np::dtype::get_builtin<float>()), is_pixels_in_log
		);
}

np::ndarray reconstruct_sirt(np::ndarray pixels,
				       np::ndarray angles,
				       double rotation_centre, int resolution,
				       int niterations, int nthreads, bool is_pixels_in_log)
{
  return reconstruct_iter(pixels, angles, rotation_centre, resolution,
			  niterations, nthreads, CCPi::alg_SIRT, 0.0, np::zeros(bp::make_tuple(1), np::dtype::get_builtin<float>()), is_pixels_in_log);
}

np::ndarray reconstruct_mlem(np::ndarray pixels,
				       np::ndarray angles,
				       double rotation_centre, int resolution,
				       int niterations, int nthreads, bool is_pixels_in_log)
{
  return reconstruct_iter(pixels, angles, rotation_centre, resolution,
			  niterations, nthreads, CCPi::alg_MLEM, 0.0, np::zeros(bp::make_tuple(1), np::dtype::get_builtin<float>()) , is_pixels_in_log);
}

np::ndarray
reconstruct_cgls_tikhonov(np::ndarray pixels,
			  np::ndarray angles,
			  double rotation_centre, int resolution,
			  int niterations, int nthreads, double regularize,
			  np::ndarray norm_r, bool is_pixels_in_log)
{
  return reconstruct_iter(pixels, angles, rotation_centre, resolution,
			  niterations, nthreads, CCPi::alg_CGLS_Tikhonov,
			  regularize, norm_r, is_pixels_in_log);
}

np::ndarray
reconstruct_cgls_tvreg(np::ndarray pixels,
		       np::ndarray angles,
		       double rotation_centre, int resolution,
		       int niterations, int nthreads,
		       double regularize, np::ndarray norm_r, bool is_pixels_in_log)
{
  return reconstruct_iter(pixels, angles, rotation_centre, resolution,
			  niterations, nthreads, CCPi::alg_CGLS_TVreg,
			  regularize, norm_r, is_pixels_in_log);
}

np::ndarray reconstruct_cgls2(np::ndarray pixels,
					np::ndarray angles,
					double rotation_centre, int resolution,
					int niterations, int nthreads,
					np::ndarray norm_r, bool is_pixels_in_log)
{
  return reconstruct_iter(pixels, angles, rotation_centre, resolution,
			  niterations, nthreads, CCPi::alg_CGLS, 0.0,
			  norm_r, is_pixels_in_log);
}

void reconstruct_tvreg()
{
}
