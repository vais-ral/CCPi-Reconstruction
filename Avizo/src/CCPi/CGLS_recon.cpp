/*
 *  Template of a compute module
 */

#include <QApplication>

#include <hxcore/HxMessage.h>
#include <hxcore/HxWorkArea.h>
#include <hxfield/HxUniformScalarField3.h>

#include "CGLS_recon.h"

#include "base_types.hpp"
#include "mpi.hpp"
#include "utils.hpp"
#include "instruments.hpp"
#include "algorithms.hpp"
#include "results.hpp"
#include "voxels.hpp"

static HxMessage *messages = 0;
static HxWorkArea *progress = theWorkArea;
static bool do_progress = false;

HX_INIT_CLASS(CGLS_recon,HxCompModule)

CGLS_recon::CGLS_recon() :
  HxCompModule(HxUniformScalarField3::getClassTypeId()),
  portAction(this, "action", QApplication::translate("CGLS_recon", "Action")),
  iterations(this, "number of iterations",
	     QApplication::translate("CGLS_recon", "Iterations")),
  resolution(this, "resolution",
	     QApplication::translate("CGLS_recon", "Pixels per Voxel")),
  beam_harden(this, "beam harden",
	      QApplication::translate("CGLS_recon", "Beam Hardening"))
{
  portAction.setLabel(0,"DoIt");
  iterations.setMinMax(5, 30);
  iterations.setValue(20);
  resolution.setMinMax(1, 8);
  resolution.setValue(1);
  beam_harden.setValue(true);
}

CGLS_recon::~CGLS_recon()
{
}

void CGLS_recon::compute()
{
  if (portAction.wasHit()) {
    do_progress = false;
    messages = theMsg;
    progress = theWorkArea;
    HxUniformScalarField3* field = (HxUniformScalarField3*) portData.source();
    // Check whether the input port is connected
    if (field == 0)
      return;
    // Todo - check that its a float field?
    run_cgls();
    // end progress area
    if (do_progress) {
      theWorkArea->stopWorking();
      theWorkArea->undivide();
    }
  }
}

void CGLS_recon::run_cgls()
{
  HxUniformScalarField3* field = (HxUniformScalarField3*) portData.source();
  const int *fdims = field->lattice.dims();
  boost::multi_array_ref<float, 3>
    pixels((float *)field->lattice.dataPtr(),
	   boost::extents[fdims[0]][fdims[1]][fdims[2]],
	   boost::fortran_storage_order());
  // Todo - fix
  boost::multi_array<float, 1> angles(boost::extents[1]);
  boost::multi_array<float, 1> h_offsets(boost::extents[1]);
  boost::multi_array<float, 1> v_offsets(boost::extents[1]);
  real source_x = -1.0;
  real detector_x = 10.0;
  real h_size = 0.1;
  real v_size = 0.1;
  // end Todo
  int pixels_per_voxel = resolution.getValue();
  int niterations = iterations.getValue();
  bool beam_hardening = beam_harden.getValue();
  CCPi::instrument *instrument = new CCPi::Nikon_XTek();
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
    // fortran order for Avizo
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
  messages->error(QApplication::translate("CGLS", message.c_str()));
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
  progress->startWorkingNoStop(QApplication::translate("CGLS", label));
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
