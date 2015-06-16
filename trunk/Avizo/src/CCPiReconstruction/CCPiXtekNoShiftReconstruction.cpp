/*
*  Template of a compute module
*/

#include <QApplication>

#include <hxcore/HxMessage.h>
#include <hxcore/HxWorkArea.h>
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
	beam_harden(this, "beam harden",
	QApplication::translate("CCPiXtekNoShiftReconstruction", "Beam Hardening"))
{
	portAction.setLabel(0,"DoIt");
	algorithm.setNum(3);
	algorithm.setLabel(0, "CGLS");
	algorithm.setLabel(1, "SIRT");
	algorithm.setLabel(2, "MLEM");
	algorithm.setValue(0);
	iterations.setMinMax(5, 30);
	iterations.setValue(20);
	beam_harden.setValue(false);
}

CCPiXtekNoShiftReconstruction::~CCPiXtekNoShiftReconstruction()
{
}

void CCPiXtekNoShiftReconstruction::compute()
{
	// Check whether the action port button was clicked
	if (!portAction.wasHit()) return;
	ccpi_recon::do_progress = false;
	ccpi_recon::messages = theMsg;
	ccpi_recon::progress = theWorkArea;
	HxUniformScalarField3* field = (HxUniformScalarField3*) portData.source();
	// Check whether the input port is connected
	if (field == 0)
	{
		theMsg->stream()<<"Is not connected to uniform scalar field source"<<std::endl;
		return;
	}
	double *angles = new double[field->lattice.dims()[2]];
	int result = field->parameters.findReal("Angles", field->lattice.dims()[2], angles);
	//Check the Rotation Angle
	if (result!=1)
	{
		theMsg->stream()<<"Angles parameter not found in input data"<<result<<std::endl;
		return;
	}
	// Todo - check of optional shifts
	int *pixel_size = new int[2];
	result = field->parameters.findNum("DetectorPixelSize", 2,pixel_size);
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
	result = field->parameters.findReal("SourceToDetector", srcToObject);
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
	HxUniformScalarField3 *field = (HxUniformScalarField3 *) portData.source();
	const int *fdims = field->lattice.dims();
	boost::multi_array_ref<float, 3>
		pixels((float *)field->lattice.dataPtr(),
		boost::extents[fdims[2]][fdims[1]][fdims[0]],
		boost::fortran_storage_order());
	float *angles_float = new float[field->lattice.dims()[2]];
	field->parameters.findReal("Angles", field->lattice.dims()[2], angles_float);
	// Todo - check of optional shifts
	int *pixel_size = new int[2];
	field->parameters.findNum("DetectorPixelSize", 2,pixel_size);
	double srcToObject;
	field->parameters.findReal("SourceToObject", srcToObject);
	double srcToDetector;
	field->parameters.findReal("SourceToDetector", srcToObject);
	int resolution;
	field->parameters.findNum("Resolution", resolution);
	boost::multi_array_ref<float, 1> angles((float *)angles_float,
		boost::extents[fdims[2]]);

	// Todo - offsets
	boost::multi_array<float, 1> h_offsets(boost::extents[1]);
	boost::multi_array<float, 1> v_offsets(boost::extents[1]);

	real h_size = pixel_size[0];
	real v_size = pixel_size[1];

	real source_x = - srcToObject;
	real detector_x = srcToDetector + source_x;

	int pixels_per_voxel = resolution;
	int niterations = iterations.getValue();
	bool beam_hardening = beam_harden.getValue();

	CCPi::instrument *instrument = 0;
	// Todo - need a type that supports this and will probably not be able
	// to use the find_centre algorithm with it.
	// if (optional shifts set)
	// instrument = ?;
	// else
	instrument = new CCPi::Nikon_XTek();

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
	}

	machine::initialise(0);
	// instrument setup from pixels/angles will probably copy
	real vox_origin[3];
	real vox_size[3];
	voxel_data *voxels = reconstruct(instrument, recon_algorithm, pixels, angles,
		h_offsets, v_offsets, pixels_per_voxel,
		source_x, detector_x, h_size, v_size,
		beam_hardening, vox_origin, vox_size);
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
		HxUniformCoord3 *coords =(HxUniformCoord3 *) output->lattice.coords();
		float *bx = coords->bbox();
		bx[0] = vox_origin[2];
		bx[1] = vox_origin[2] + float(dims[0]) * vox_size[2];
		bx[2] = vox_origin[1];
		bx[3] = vox_origin[1] + float(dims[0]) * vox_size[1];
		bx[4] = vox_origin[0];
		bx[5] = vox_origin[0] + float(dims[2]) * vox_size[0];
		// publish reconstruction
		setResult(output); 
	}
}
