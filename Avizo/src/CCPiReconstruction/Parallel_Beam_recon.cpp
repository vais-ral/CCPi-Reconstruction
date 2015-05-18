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

static HxMessage *messages = 0;
static HxWorkArea *progress = theWorkArea;
static bool do_progress = false;

HX_INIT_CLASS(Parallel_Beam_recon,HxCompModule)

Parallel_Beam_recon::Parallel_Beam_recon() :
  HxCompModule(HxUniformScalarField3::getClassTypeId()),
  portAction(this, "action", QApplication::translate("Parallel_Beam_recon", "Action")),
  rotationAngle(this,"rotation angle",QApplication::translate("Parallel_Beam_recon", "Rotation Angle")),
  pixelSize(this,"pixel size",QApplication::translate("Parallel_Beam_recon", "Pixel Size(x,y)")),
  imageKey(this,"image key",QApplication::translate("Parallel_Beam_recon", "Image Key")),
  iterations(this, "number of iterations",
	     QApplication::translate("Parallel_Beam_recon", "Iterations")),
  resolution(this, "resolution",
	     QApplication::translate("Parallel_Beam_recon", "Pixels per Voxel")),
  beam_harden(this, "beam harden",
	      QApplication::translate("Parallel_Beam_recon", "Beam Hardening"))
{
  rotationAngle.addType(HxUniformScalarField3::getClassTypeId());
  pixelSize.addType(HxUniformScalarField3::getClassTypeId());
  imageKey.addType(HxUniformScalarField3::getClassTypeId());

  portAction.setLabel(0,"DoIt");
  iterations.setMinMax(5, 30);
  iterations.setValue(20);
  resolution.setMinMax(1, 8);
  resolution.setValue(1);
  beam_harden.setValue(true);
}

Parallel_Beam_recon::~Parallel_Beam_recon()
{
}

void Parallel_Beam_recon::compute()
{
  if (portAction.wasHit()) {
    do_progress = false;
    messages = theMsg;
    progress = theWorkArea;
    HxUniformScalarField3* field = (HxUniformScalarField3*) portData.source();
    // Check whether the input port is connected
    if (field == 0)
      return;
	HxUniformScalarField3* rot_angle = (HxUniformScalarField3*) rotationAngle.source();
	//Check the Rotation Angle
	if (rot_angle == 0 || rot_angle->lattice.dims()[0]==1)
		return;
	HxUniformScalarField3* image_key = (HxUniformScalarField3*) imageKey.source();
	//Check the image key
	if (image_key == 0 || image_key->lattice.dims()[0]==1)
		return;


	
    // Todo - check that its a float field?
    run_reconstruction();
    // end progress area
    if (do_progress) {
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
  messages->error(QApplication::translate("Parallel_Beam",std::to_string((long double)fdims[2]).c_str()));// std::to_string((long double)fdims[2])));
  float* tmp_angles = (float*)((HxUniformScalarField3*) rotationAngle.source())->lattice.dataPtr();
  float* image_key = (float*)((HxUniformScalarField3*) imageKey.source())->lattice.dataPtr();
  // Todo - fix
  boost::multi_array<float, 1> angles(boost::extents[fdims[2]]);
  boost::multi_array<float, 1> h_offsets(boost::extents[1]);
  boost::multi_array<float, 1> v_offsets(boost::extents[1]);

  for(int i=0;i<fdims[2];i++)
	  angles[i] = tmp_angles[i];
  real source_x = -1.0;
  real detector_x = 10.0;
  real h_size = 0.1;
  real v_size = 0.1;
  // end Todo
  int pixels_per_voxel = resolution.getValue();
  int niterations = iterations.getValue();
  bool beam_hardening = beam_harden.getValue();
  CCPi::instrument *instrument = new CCPi::Diamond();
  CCPi::reconstruction_alg *algorithm = new CCPi::cgls_3d(niterations);
  //if (blocking_factor > 0 and instrument->supports_blocks())
  //  recon_algorithm = new CCPi::cgls_2d(niterations, pixels_per_voxel);
  machine::initialise(0);
  // instrument setup from pixels/angles will probably copy
  voxel_data *voxels = reconstruct(instrument, algorithm, pixels, angles,
				   h_offsets, v_offsets, pixels_per_voxel,
				   source_x, detector_x, h_size, v_size,
				   beam_hardening);
  machine::exit();
  delete algorithm;
  delete instrument;
  if (voxels != 0) {
    int dims[3];
    
    dims[0] = voxels->shape()[2];
    dims[1] = voxels->shape()[1];
    dims[2] = voxels->shape()[0];
    HxUniformScalarField3* output =
      new HxUniformScalarField3(dims, field->primType());
    for (int i = 0; i < dims[0]; i++)
      for (int j = 0; j < dims[1]; j++)
	for (int k = 0; k < dims[2]; k++)
	  output->set(i, j, k, (*voxels)[i][j][k]);
    delete voxels;
    // Todo - some sort of setBoundingBox for output? based on origin etc
    setResult(output); 
  }
}

void report_error(const std::string message)
{
  messages->error(QApplication::translate("Parallel_Beam", message.c_str()));
}

void report_error(const std::string message, const std::string arg)
{
  std::string sum = message + arg;
  report_error(sum);
}

void report_error(const std::string message, const std::string arg1,
		  const std::string arg2)
{
  std::string sum = message + arg1;
  sum += arg2;
  report_error(sum);
}

void initialise_progress(const int length, const char label[])
{
  if (do_progress) {
    progress->stopWorking();
    progress->undivide();
  }
  // Turn into busy state, don't activate the Stop button.
  progress->subdivide(length);
  progress->startWorkingNoStop(QApplication::translate("Parallel_Beam", label));
  do_progress = true;
}

void update_progress(const int value)
{
  //progress->progressStep();
  progress->setProgressValue(value);
}

void add_output(const std::string str)
{
  messages->printf("%s", str.c_str());
}

void add_output(const char c)
{
  messages->printf("%c", c);
}

void add_output(const int i)
{
  messages->printf("%1d", i);
}

void add_output(const sl_int i)
{
  messages->printf("%1ld", i);
}

void add_output(const int i, const int w, const bool fill)
{
  if (fill)
    messages->printf("%0*d", w, i);
  else
    messages->printf("%*d", w, i);
}

void add_output(const real r)
{
  messages->printf("%f", r);
}

void send_output()
{
  messages->printf("\n");
}
