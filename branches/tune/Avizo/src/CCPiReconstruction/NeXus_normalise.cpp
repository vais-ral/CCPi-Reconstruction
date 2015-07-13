/*
 *  Template of a compute module
 */

// This is basically a test module for the parallel beam reconstruction since
// we don't have all the filters etc that we would use in SAVU.

#include <QApplication>

#include <hxcore/HxMessage.h>
#include <hxcore/HxWorkArea.h>
#include <hxfield/HxUniformScalarField3.h>

#include "NeXus_normalise.h"
#include "base_types.hpp"

HX_INIT_CLASS(NeXus_normalise,HxCompModule)

NeXus_normalise::NeXus_normalise()
: HxCompModule(HxUniformScalarField3::getClassTypeId()),
  portAction(this,"action",QApplication::translate("NeXus_normalise", "Action")),
  rotationAngle(this,"rotation angle",QApplication::translate("NeXus_normalise", "Rotation Angle")),
  imageKey(this,"image key",QApplication::translate("NeXus_normalise", "Image Key")),
  pixelSize(this,"pixel size",QApplication::translate("NeXus_normalise", "Pixel Size(x,y)"))
{
  rotationAngle.addType(HxUniformScalarField3::getClassTypeId());
  imageKey.addType(HxUniformScalarField3::getClassTypeId());
  pixelSize.addType(HxUniformScalarField3::getClassTypeId());
  portAction.setLabel(0,"DoIt");
}

NeXus_normalise::~NeXus_normalise()
{
}

void NeXus_normalise::compute()
{
  if (portAction.wasHit()) {
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
    //Check the image key
    HxUniformScalarField3* image_key =
      (HxUniformScalarField3*) imageKey.source();
    if (image_key == 0 or image_key->lattice.dims()[0] == 1)
      return;
    // Is this the right order?
    if (image_key->lattice.dims()[0] != field->lattice.dims()[2]) {
      theMsg->error(QApplication::translate("NeXus_normalise",
					    "Image key/data dim mismatch"));
      return;
    }
    // check types - are these right?
    if (image_key->lattice.primType() != McPrimType::mc_int32) {
      theMsg->error(QApplication::translate("NeXus_normalise",
					    "Incorrect type for image keys"));
      return;
    }
    if (field->lattice.primType() != McPrimType::mc_int16) {
      theMsg->error(QApplication::translate("NeXus_normalise",
				  "Normalisation of data type not supported"));
      return;
    } else
      normalise();
  }
}

void NeXus_normalise::normalise()
{
  HxUniformScalarField3* field = (HxUniformScalarField3*) portData.source();
  const int *fdims = field->lattice.dims();
  boost::multi_array_ref<short, 3>
    data((short *)field->lattice.dataPtr(),
	 boost::extents[fdims[0]][fdims[1]][fdims[2]],
	 boost::fortran_storage_order());
  HxUniformScalarField3* image_key = (HxUniformScalarField3*) imageKey.source();
  int *keys = (int *)image_key->lattice.dataPtr();
  int nkeys = image_key->lattice.dims()[0];
  HxUniformScalarField3* old_angles =
    (HxUniformScalarField3*) rotationAngle.source();
  // Todo - is this the correct type
  float *angles = (float *)old_angles->lattice.dataPtr();

  // find number of projections
  int nprojections = 0;
  for (int i = 0; i < nkeys; i++)
    if (keys[i] == 0)
      nprojections++;
  // create field of correct size
  int dims[3];
  dims[0] = fdims[0];
  dims[1] = fdims[1];
  dims[2] = nprojections;
  HxUniformScalarField3* output =
    new HxUniformScalarField3(dims, McPrimType::mc_float);
  boost::multi_array_ref<float, 3>
    pixels((float *)output->lattice.dataPtr(),
	   boost::extents[dims[0]][dims[1]][dims[2]],
	   boost::fortran_storage_order());
  // update angles
  int adims[3];
  adims[0] = nprojections;
  adims[1] = 1;
  adims[2] = 2;
  HxUniformScalarField3* new_angles =
    new HxUniformScalarField3(adims, McPrimType::mc_float);
  float *aptr = (float *)new_angles->lattice.dataPtr();
  int proj_index = 0;
  int n_ibright = 0;
  int n_fbright = 0;
  int n_idark = 0;
  int n_fdark = 0;
  // boost initialises these so should be ok
  boost::multi_array<float, 2> i_bright(boost::extents[fdims[0]][fdims[1]],
					boost::fortran_storage_order());
  boost::multi_array<float, 2> f_bright(boost::extents[fdims[0]][fdims[1]],
					boost::fortran_storage_order());
  boost::multi_array<float, 2> i_dark(boost::extents[fdims[0]][fdims[1]],
				      boost::fortran_storage_order());
  boost::multi_array<float, 2> f_dark(boost::extents[fdims[0]][fdims[1]],
				      boost::fortran_storage_order());
  for (int k = 0; k < fdims[2]; k++) {
    if (keys[k] == 0) {
      for (int j = 0; j < fdims[1]; j++) {
	for (int i = 0; i < fdims[0]; i++) {
	  pixels[i][j][proj_index] = (float)data[i][j][k];
	}
      }
      aptr[proj_index] = angles[k];
      proj_index++;
    } else if (keys[k] == 1) {
      if (k > fdims[2] / 2) {
	for (int j = 0; j < fdims[1]; j++) {
	  for (int i = 0; i < fdims[0]; i++) {
	    f_bright[i][j] += (float)data[i][j][k];
	  }
	}
	n_fbright++;
      } else {
	for (int j = 0; j < fdims[1]; j++) {
	  for (int i = 0; i < fdims[0]; i++) {
	    i_bright[i][j] += (float)data[i][j][k];
	  }
	}
	n_ibright++;
      }
    } else if (keys[k] == 2) {
      if (k > fdims[2] / 2) {
	for (int j = 0; j < fdims[1]; j++) {
	  for (int i = 0; i < fdims[0]; i++) {
	    f_dark[i][j] += (float)data[i][j][k];
	  }
	}
	n_fdark++;
      } else {
	for (int j = 0; j < fdims[1]; j++) {
	  for (int i = 0; i < fdims[0]; i++) {
	    i_dark[i][j] += (float)data[i][j][k];
	  }
	}
	n_idark++;
      }
    }
  }
  // we're averaging the flats/darks
  if (n_idark > 0) {
    for (int j = 0; j < fdims[1]; j++)
      for (int i = 0; i < fdims[0]; i++)
	i_dark[i][j] /= float(n_idark);
  }
  if (n_fdark > 0) {
    for (int j = 0; j < fdims[1]; j++)
      for (int i = 0; i < fdims[0]; i++)
	f_dark[i][j] /= float(n_fdark);
  }
  if (n_ibright > 0) {
    for (int j = 0; j < fdims[1]; j++)
      for (int i = 0; i < fdims[0]; i++)
	i_bright[i][j] /= float(n_ibright);
  }
  if (n_fbright > 0) {
    for (int j = 0; j < fdims[1]; j++)
      for (int i = 0; i < fdims[0]; i++)
	f_bright[i][j] /= float(n_fbright);
  }
  if (n_ibright == 0) {
    if (n_fbright == 0) {
      //add_output("Warning no flats found");
      //send_output();
      // find max value and use that?
      float max_pixel = 0.0;
      for (int k = 0; k < nprojections; k++) {
	for (int j = 0; j < fdims[1]; j++) {
	  for (int i = 0; i < fdims[0]; i++) {
	    if (pixels[i][j][k] > max_pixel)
	      max_pixel = pixels[i][j][k];
	  }
	}
      }
      for (int j = 0; j < fdims[1]; j++)
	for (int i = 0; i < fdims[0]; i++)
	  i_bright[i][j] = max_pixel;
      for (int j = 0; j < fdims[1]; j++)
	for (int i = 0; i < fdims[0]; i++)
	  f_bright[i][j] = max_pixel;
    } else {
      for (int j = 0; j < fdims[1]; j++)
	for (int i = 0; i < fdims[0]; i++)
	  i_bright[i][j] = f_bright[i][j];
    }
  } else if (n_fbright == 0) {
    for (int j = 0; j < fdims[1]; j++)
      for (int i = 0; i < fdims[0]; i++)
	f_bright[i][j] = i_bright[i][j];
  }
  // we interpolate the flats/darks from each end to the angle.
  boost::multi_array<float, 2> dark(boost::extents[fdims[0]][fdims[1]],
				    boost::fortran_storage_order());
  boost::multi_array<float, 2> bright(boost::extents[fdims[0]][fdims[1]],
				      boost::fortran_storage_order());
  for (int k = 0; k < nprojections; k++) {
    // Based on fbp code, interpolate bright/dark
    real w = (aptr[k] - aptr[0]) / (aptr[nprojections - 1] - aptr[0]);
    for (int j = 0; j < fdims[1]; j++)
      for (int i = 0; i < fdims[0]; i++)
	dark[i][j] = i_dark[i][j] * (1.0 - w) + f_dark[i][j] * w;
    for (int j = 0; j < fdims[1]; j++)
      for (int i = 0; i < fdims[0]; i++)
	bright[i][j] = i_bright[i][j] * (1.0 - w) + f_bright[i][j] * w;
    // subtract dark from data/bright
    // and clamp min data/bright value to 0.1
    for (int j = 0; j < fdims[1]; j++) {
      for (int i = 0; i < fdims[0]; i++) {
	bright[i][j] -= dark[i][j];
	if (bright[i][j] < 0.1)
	  bright[i][j] = 0.1;
      }
    }
    for (int j = 0; j < fdims[1]; j++) {
      for (int i = 0; i < fdims[0]; i++) {
	pixels[i][j][k] -= dark[i][j];
	if (pixels[i][j][k] < 0.1)
	  pixels[i][j][k] = 0.1;
      }
    }
    // scale each data pixel by bright pixel
    for (int j = 0; j < fdims[1]; j++)
      for (int i = 0; i < fdims[0]; i++)
	pixels[i][j][k] /= bright[i][j];
    // clamp to 1.0
    for (int j = 0; j < fdims[1]; j++) {
      for (int i = 0; i < fdims[0]; i++) {
	if (pixels[i][j][k] > 1.0)
	  pixels[i][j][k] = 1.0;
      }
    }
  }
  // register result
  HxData::registerData(output, "projections");
  HxData::registerData(new_angles, "rotation angles");
  // pass through pixels - is this the correct way to do it?
  HxData::registerData((HxData *)pixelSize.source(), "pixels (x,y)");
}
