/*
 *  Template of a compute module
 */

#include <QApplication>

#include <hxcore/HxMessage.h>
#include <hxcore/HxWorkArea.h>
#include <hxfield/HxUniformScalarField3.h>

#include "Parallel_Beam_recon.h"

#include "base_types.hpp"
#include "mpi.hpp"
#include "utils.hpp"
#include "instruments.hpp"
#include "algorithms.hpp"
#include "cgls.hpp"
#include "sirt.hpp"
#include "results.hpp"
#include "voxels.hpp"
#include "ui_messages.hpp"

HX_INIT_CLASS(Parallel_Beam_recon,HxCompModule)

Parallel_Beam_recon::Parallel_Beam_recon() :
  HxCompModule(HxUniformScalarField3::getClassTypeId()),
  portAction(this, "action", QApplication::translate("Parallel_Beam_recon",
						     "Action")),
  rotationAngle(this,"rotation angle",
		QApplication::translate("Parallel_Beam_recon",
					"Rotation Angle")),
  pixelSize(this,"pixel size",
	    QApplication::translate("Parallel_Beam_recon", "Pixel Size(x,y)")),
  algorithm(this, "algorithm",
	    QApplication::translate("Parallel_Beam_recon",
				    "Reconstruction Algorithm")),
  iterations(this, "number of iterations",
	     QApplication::translate("Parallel_Beam_recon", "Iterations")),
  resolution(this, "resolution",
	     QApplication::translate("Parallel_Beam_recon",
				     "Pixels per Voxel")),
  rotationCentre(this, "rotation centre",
		 QApplication::translate("Parallel_Beam_recon",
					 "Rotation Centre")),
  beam_harden(this, "beam harden",
	      QApplication::translate("Parallel_Beam_recon", "Beam Hardening"))
{
  rotationAngle.addType(HxUniformScalarField3::getClassTypeId());
  pixelSize.addType(HxUniformScalarField3::getClassTypeId());

  portAction.setLabel(0,"DoIt");
  algorithm.setNum(2);
  algorithm.setLabel(0, "CGLS");
  algorithm.setLabel(0, "SIRT");
  algorithm.setValue(0);
  iterations.setMinMax(5, 30);
  iterations.setValue(20);
  resolution.setMinMax(1, 8);
  resolution.setValue(1);
  rotationCentre.setNum(1);
  rotationCentre.setValue(0, 0.0);
  beam_harden.setValue(false);
}

Parallel_Beam_recon::~Parallel_Beam_recon()
{
}

void Parallel_Beam_recon::compute()
{
  if (portAction.wasHit()) {
    ccpi_recon::do_progress = false;
    ccpi_recon::messages = theMsg;
    ccpi_recon::progress = theWorkArea;
    HxUniformScalarField3* field = (HxUniformScalarField3*) portData.source();
    // Check whether the input port is connected
    if (field == 0)
      return;
    HxUniformScalarField3* rot_angle =
      (HxUniformScalarField3*) rotationAngle.source();
    //Check the Rotation Angle
    if (rot_angle == 0 or rot_angle->lattice.dims()[0] == 1)
      return;
    //Check the Pixels
    HxUniformScalarField3* pixel_size =
      (HxUniformScalarField3*) pixelSize.source();
    if (pixel_size == 0 or pixel_size->lattice.dims()[0] != 2)
      return;
	
    // Todo - check that its a float field?
    run_reconstruction();
    // end progress area
    if (ccpi_recon::do_progress) {
      theWorkArea->stopWorking();
      theWorkArea->undivide();
    }
  }
}

void Parallel_Beam_recon::run_reconstruction()
{
  HxUniformScalarField3* field = (HxUniformScalarField3*) portData.source();
  const int *fdims = field->lattice.dims();
  boost::multi_array_ref<float, 3>
    pixels((float *)field->lattice.dataPtr(),
	   boost::extents[fdims[0]][fdims[1]][fdims[2]],
	   boost::fortran_storage_order());
  //messages->error(QApplication::translate("Parallel_Beam",std::to_string((long double)fdims[2]).c_str()));// std::to_string((long double)fdims[2])));
  HxUniformScalarField3* rot_angle =
    (HxUniformScalarField3*) rotationAngle.source();
  boost::multi_array_ref<float, 1> angles((float *)rot_angle->lattice.dataPtr(),
					  boost::extents[fdims[2]]);

  HxUniformScalarField3* pix = (HxUniformScalarField3*) pixelSize.source();
  float *pixel_ptr = (float *)pix->lattice.dataPtr();
  real h_size = pixel_ptr[0];
  real v_size = pixel_ptr[1];

  real rotation_centre = rotationCentre.getValue();
  int pixels_per_voxel = resolution.getValue();
  int niterations = iterations.getValue();
  bool beam_hardening = beam_harden.getValue();
  CCPi::instrument *instrument = new CCPi::Diamond();

  CCPi::reconstruction_alg *recon_algorithm = 0;
  switch (algorithm.getValue()) {
  case 0:
    //if (blocking_factor > 0 and instrument->supports_blocks())
    //  recon_algorithm = new CCPi::cgls_2d(niterations, pixels_per_voxel);
    recon_algorithm = new CCPi::cgls_3d(niterations);
    break;
  case 1:
    recon_algorithm = new CCPi::sirt(niterations);
    break;
  }

  machine::initialise(0);
  // instrument setup from pixels/angles will probably copy
  voxel_data *voxels = reconstruct(instrument, recon_algorithm, pixels, angles,
				   rotation_centre, pixels_per_voxel, 0,
				   beam_hardening);
  machine::exit();
  delete recon_algorithm;
  delete instrument;
  if (voxels != 0) {
    int dims[3];
    
    dims[0] = voxels->shape()[2];
    dims[1] = voxels->shape()[1];
    dims[2] = voxels->shape()[0];
    HxUniformScalarField3* output =
      new HxUniformScalarField3(dims, McPrimType::mc_float);
    for (int i = 0; i < dims[0]; i++)
      for (int j = 0; j < dims[1]; j++)
	for (int k = 0; k < dims[2]; k++)
	  output->set(i, j, k, (*voxels)[i][j][k]);
    delete voxels;
    // set bounding box for coords based on pixel sizes - is this right?
    HxUniformCoord3 *coords =(HxUniformCoord3 *) output->lattice.coords();
    float *bx = coords->bbox();
    // z extent - (dims[0] / 2) * (pixels_per_voxel * v_size)
    bx[1] = float(dims[0] / 2) * pixels_per_voxel * v_size;
    bx[0] = - bx[1];
    // x/y extent, should be the same
    bx[3] = float(dims[1] / 2) * pixels_per_voxel * h_size;
    bx[2] = - bx[1];
    bx[5] = float(dims[2] / 2) * pixels_per_voxel * h_size;
    bx[4] = - bx[1];
    // publish reconstruction
    setResult(output); 
  }
}
