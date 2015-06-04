/*
 *  Template of a compute module
 */

#include <QApplication>

#include <hxcore/HxMessage.h>
#include <hxcore/HxWorkArea.h>
#include <hxfield/HxUniformScalarField3.h>

#include "XTek_recon.h"
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

HX_INIT_CLASS(XTek_recon,HxCompModule)

XTek_recon::XTek_recon() :
  HxCompModule(HxUniformScalarField3::getClassTypeId()),
  portAction(this,"action",QApplication::translate("XTek_recon", "Action")),
  algorithm(this, "algorithm",
	    QApplication::translate("XTek_recon",
				    "Reconstruction Algorithm")),
  iterations(this, "number of iterations",
	     QApplication::translate("XTek_recon", "Iterations")),
  resolution(this, "resolution",
	     QApplication::translate("XTek_recon", "Pixels per Voxel")),
  beam_harden(this, "beam harden",
	      QApplication::translate("XTek_recon", "Beam Hardening"))
{
  portAction.setLabel(0,"DoIt");
  algorithm.setNum(3);
  algorithm.setLabel(0, "CGLS");
  algorithm.setLabel(0, "SIRT");
  algorithm.setLabel(0, "MLEM");
  algorithm.setValue(0);
  iterations.setMinMax(5, 30);
  iterations.setValue(20);
  resolution.setMinMax(1, 8);
  resolution.setValue(1);
  beam_harden.setValue(false);
}

XTek_recon::~XTek_recon()
{
}

void XTek_recon::compute()
{
  if (portAction.wasHit()) {
    HxFileDialog *dialog = HxFileDialog::_getTheFileDialog();
    QString filter = dialog->getFileNameFilter(QString::fromAscii("XTek"),
					       QString::fromAscii("xketct"));
    QString *selected_filter;
    QString filename;
    QString thisfilter;
    dialog->getFileNameAndFilter(QString::fromAscii("XTek Data"),
				 QString::fromAscii(""), filter,
				 selected_filter, filename, thisfilter);
    ccpi_recon::do_progress = false;
    ccpi_recon::messages = theMsg;
    ccpi_recon::progress = theWorkArea;
    const QByteArray asc = filename.toAscii();
    std::string file(asc.constData(), asc.length());
    run_reconstruction(file);
    // end progress area
    if (ccpi_recon::do_progress) {
      theWorkArea->stopWorking();
      theWorkArea->undivide();
    }
  }
}


void XTek_recon::run_reconstruction(const std::string filename)
{
  // Todo - need filename from browser
  int pixels_per_voxel = resolution.getValue();
  int niterations = iterations.getValue();
  bool beam_hardening = beam_harden.getValue();

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
  }

  machine::initialise(0);
  // instrument setup from pixels/angles will probably copy
  real vox_origin[3];
  real vox_size[3];
  voxel_data *voxels = reconstruct(instrument, recon_algorithm, filename, "",
				   vox_origin, vox_size, -1.0,
				   pixels_per_voxel, 0, beam_hardening,
				   CCPi::no_output, false, false);
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
