/*
*  Template of a compute module
*/

#include <QApplication>

#include <hxcore/HxMessage.h>
#include <hxcore/internal/HxWorkArea.h>
#include <hxfield/HxUniformScalarField3.h>

#include "CCPiXtekNoShiftReconstruction.h"
#include "base_types.hpp"
#include "mpi.hpp"
#include "utils.hpp"
#include "instruments.hpp"
#include "algorithms.hpp"
#include "cgls.hpp"
#include "sirt.hpp"
#include "mlem.hpp"
#include "results.hpp"
#include "voxels.hpp"
#include "ui_messages.hpp"
#include <boost/multi_array.hpp>

HX_INIT_CLASS(CCPiXtekNoShiftReconstruction, HxCompModule)

	CCPiXtekNoShiftReconstruction::CCPiXtekNoShiftReconstruction() :
HxCompModule(HxUniformScalarField3::getClassTypeId()),
	portAction(this,"action",QApplication::translate("CCPiXtekNoShiftReconstruction",
	"Action")),
	algorithm(this, "algorithm",
	QApplication::translate("CCPiXtekNoShiftReconstruction",
	"Reconstruction Algorithm")),
	iterations(this, "number of iterations",
	QApplication::translate("CCPiXtekNoShiftReconstruction", "Iterations")),
	resolution(this, "resolution", QApplication::translate("CCPiXtekNoShiftReconstruction", "Pixels per Voxel")),
	beam_harden(this, "beam harden",
	QApplication::translate("CCPiXtekNoShiftReconstruction", "Beam Hardening")),
	regularize(this, "regularize", QApplication::translate("CCPiXtekNoShiftReconstruction", "Regularization parameter"))
{
	portAction.setLabel(0,"DoIt");
	algorithm.setNum(5);
	algorithm.setLabel(0, "CGLS");
	algorithm.setLabel(1, "SIRT");
	algorithm.setLabel(2, "MLEM");
	algorithm.setLabel(3, "CGLS/Tikhonov");
	algorithm.setLabel(4, "CGLS/TV reg");
	algorithm.setValue(0);
	iterations.setMinMax(5, 30);
	iterations.setValue(20);
	resolution.setMinMax(1, 8);
	resolution.setValue(1);
	beam_harden.setValue(false);
	regularize.setNum(1);
	regularize.setMinMax(0.00001, 0.2);
	regularize.setValue(0.01);
	regularize.hide();
}

CCPiXtekNoShiftReconstruction::~CCPiXtekNoShiftReconstruction()
{
}

void CCPiXtekNoShiftReconstruction::compute()
{
	// Check whether the action port button was clicked
  if (!portAction.wasHit()) {
    if (algorithm.isNew()) {
      if (algorithm.getValue() < 3)
	regularize.hide();
      else
	regularize.show();
    }
    return;
  }
	ccpi_recon::do_progress = false;
	ccpi_recon::messages = theMsg;
	ccpi_recon::progress = theWorkArea;
	HxUniformScalarField3* field = (HxUniformScalarField3*) portData.getSource();
	// Check whether the input port is connected
	if (field == 0)
	{
		theMsg->stream()<<"Is not connected to uniform scalar field source"<<std::endl;
		return;
	}
	float *angles = new float[field->lattice().getDims()[2]];
	int result = field->parameters.findReal("Angles", field->lattice().getDims()[2], angles);
	//Check the Rotation Angle
	if (result!=1)
	{
		theMsg->stream()<<"Angles parameter not found in input data"<<result<<std::endl;
		return;
	}
	// Todo - check of optional shifts
	float *pixel_size = new float[2];
	result = field->parameters.findReal("DetectorPixelSize", 2,pixel_size);
	if(result!=1)
	{
		theMsg->stream()<<"DetectorPixelSize parameter not found in input data"<<result<<std::endl;
		return;
	}
	double srcToObject;
	result = field->parameters.findReal("SourceToObject", srcToObject);
	if(result!=1)
	{
		theMsg->stream()<<"SourceToObject parameter not found in input data"<<std::endl;
		return;
	}
	double srcToDetector;
	result = field->parameters.findReal("SourceToDetector", srcToDetector);
	if(result!=1)
	{
		theMsg->stream()<<"SourceToDetector parameter not found in input data"<<std::endl;
		return;
	}

	// Todo - check that its a float field?
	run_reconstruction();

	// end progress area
	if (ccpi_recon::do_progress) {
		theWorkArea->stopWorking();
		theWorkArea->undivide();
	}
}

void CCPiXtekNoShiftReconstruction::run_reconstruction()
{
	HxUniformScalarField3 *field = (HxUniformScalarField3 *) portData.getSource();
	McDim3l fdims = field->lattice().getDims();
	boost::multi_array_ref<float, 3>
		pixels_ref((float *)field->lattice().dataPtr(),
		boost::extents[fdims[0]][fdims[1]][fdims[2]],
		boost::fortran_storage_order());

	//Convert pixels to correct order
	pixel_3d pixels(boost::extents[fdims[2]][fdims[1]][fdims[0]]);
	for (int i = 0; i < fdims[2]; i++)
		for (int j = 0; j < fdims[1]; j++)
			for (int k = 0; k < fdims[0]; k++)
				pixels[i][j][k] = pixels_ref[k][j][i];

	float *angles_float = new float[field->lattice().getDims()[2]];
	field->parameters.findReal("Angles", field->lattice().getDims()[2], angles_float);
	float *pixel_size = new float[2];
	field->parameters.findReal("DetectorPixelSize", 2,pixel_size);
	double srcToObject;
	field->parameters.findReal("SourceToObject", srcToObject);
	double srcToDetector;
	field->parameters.findReal("SourceToDetector", srcToDetector);
	double MaskRadius;
	field->parameters.findReal("MaskRadius", MaskRadius);
	boost::multi_array_ref<float, 1> angles(angles_float,
						boost::extents[fdims[2]]);

	real h_size = pixel_size[0];
	real v_size = pixel_size[1];

	real source_x = - srcToObject;
	real detector_x = srcToDetector + source_x;

	int pixels_per_voxel = resolution.getValue();
	int niterations = iterations.getValue();
	bool beam_hardening = beam_harden.getValue();
	float reg_param = regularize.getValue();

	CCPi::instrument *instrument = new CCPi::Nikon_XTek();

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
	case 2:
		recon_algorithm = new CCPi::mlem(niterations);
		break;
	case 3:
	  recon_algorithm = new CCPi::cgls_tikhonov(niterations, reg_param);
		break;
	case 4:
	  recon_algorithm = new CCPi::cgls_tv_reg(niterations, reg_param);
		break;
	}

	machine::initialise(0);
	// instrument setup from pixels/angles will probably copy
	boost::multi_array<float, 1> h_offsets(boost::extents[1]);
	boost::multi_array<float, 1> v_offsets(boost::extents[1]);
	real vox_origin[3];
	real vox_size[3];
	voxel_data *voxels = reconstruct(instrument, recon_algorithm, pixels,
					 angles, h_offsets, v_offsets,
					 pixels_per_voxel, source_x, detector_x,
					 h_size, v_size, MaskRadius,
					 beam_hardening, vox_origin, vox_size,
					 false, false);
	theMsg->stream() << "completed the reconstruction: " << voxels<< std::endl;

	machine::exit();
	delete recon_algorithm;
	delete instrument;
	if (voxels != 0) {
		int dims[3];

		dims[0] = voxels->shape()[0];
		dims[1] = voxels->shape()[1];
		dims[2] = voxels->shape()[2];
		HxUniformScalarField3* output =
			new HxUniformScalarField3(dims, McPrimType::MC_FLOAT);
		for (int i = 0; i < dims[0]; i++)
			for (int j = 0; j < dims[1]; j++)
				for (int k = 0; k < dims[2]; k++)
					output->set(i, j, k, (*voxels)[i][j][k]);
		delete voxels;
		HxUniformCoord3 *coords =(HxUniformCoord3 *) output->lattice().coords();
		McBox3f bx = coords->getBoundingBox();
		bx[0] = vox_origin[0];
		bx[1] = vox_origin[0] + float(dims[0]) * vox_size[0];
		bx[2] = vox_origin[1];
		bx[3] = vox_origin[1] + float(dims[1]) * vox_size[1];
		bx[4] = vox_origin[2];
		bx[5] = vox_origin[2] + float(dims[2]) * vox_size[2];
		// publish reconstruction
		setResult(output); 
	}
}
