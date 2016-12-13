/*
 *  Template of a compute module
 */

// This is basically a test module for the parallel beam reconstruction since
// we don't have all the filters etc that we would use in SAVU.

#include <QApplication>

#include <hxcore/HxMessage.h>
#include <hxcore/HxWorkArea.h>
#include <hxcore/HxObjectPool.h>
#include <hxfield/HxUniformScalarField3.h>
#include <hxfield/HxMultiChannelField3.h>

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
	  HxMultiChannelField3 *mcf = new HxMultiChannelField3();
	  mcf->setLabel("TomoProjection");
	  theObjectPool->addObject(mcf);
	  if (validateAndPopulateData(mcf))
	  {
		  normalise(mcf);
		  setResult(mcf);
	  }
	  else
		  theObjectPool->removeObject(mcf);
  }
}

bool NeXus_normalise::validateAndPopulateData(HxMultiChannelField3 *mcf)
{
	HxUniformScalarField3* field = (HxUniformScalarField3*) portData.source();
	// Check whether the input port is connected
	if (field == 0)
		return false;
	HxUniformScalarField3* rot_angle =
		(HxUniformScalarField3*) rotationAngle.source();
	//Check the Rotation Angle 
	//if angle is not assigned then check whether the field has the parameter 'Angles'
	if (rot_angle == 0 or rot_angle->lattice.dims()[0] == 1)
	{	
		HxParameter* param = field->parameters.find("Angles");
		if(param==NULL)
			return false;
		else
		{
			angles = new double[param->dim()];
			param->getReal(angles);
		}
	} else {
		angles = (double*)rot_angle->lattice.dataPtr();
	}
	//Check the Pixels
	HxUniformScalarField3* pixel_size =
		(HxUniformScalarField3*) pixelSize.source();
	if (pixel_size == 0 or pixel_size->lattice.dims()[0] != 2)
	{
		HxParameter* param = field->parameters.find("PixelSize");
		if(param==NULL)
			return false;
		else
		{
			pixelSizeXY = new double[param->dim()];
			param->getReal(pixelSizeXY);
			//Register new dataset since pixel size dataset is not available
			int adims[3];
			adims[0]=2;
			adims[1]=1;
			adims[2]=1;
			HxUniformScalarField3* pixelData = new HxUniformScalarField3(adims, McPrimType::mc_float);
			float *pixelDataPtr = (float*)pixelData->lattice.dataPtr();
			pixelDataPtr[0] = pixelSizeXY[0];
			pixelDataPtr[1] = pixelSizeXY[1];
			HxData::registerData(pixelData, "pixels (x,y)");
			pixelData->portMaster.connect(mcf);
		}
	}
	//Check the image key
	HxUniformScalarField3* image_key =
		(HxUniformScalarField3*) imageKey.source();
	if (image_key == 0 or image_key->lattice.dims()[0] == 1)
	{
		HxParameter* param = field->parameters.find("ImageKey");
		if(param==NULL)
			return false;
		else
		{
			numberOfKeys = param->dim();
			imageKeyIds = new int[numberOfKeys];
			param->getNum(imageKeyIds);
		}
	} else {
		numberOfKeys = image_key->lattice.dims()[0];
		// Is this the right order?
		if (numberOfKeys != field->lattice.dims()[2]) {
			theMsg->error(QApplication::translate("NeXus_normalise",
				"Image key/data dim mismatch"));
			return false;
		}
		// check types - are these right?
		if (image_key->lattice.primType() != McPrimType::mc_int32) {
			theMsg->error(QApplication::translate("NeXus_normalise",
				"Incorrect type for image keys"));
			return false;
		}
		imageKeyIds = (int *)image_key->lattice.dataPtr();
	}
	return true;
}

void NeXus_normalise::normalise(HxMultiChannelField3 *mcf)
{
  HxUniformScalarField3* field = (HxUniformScalarField3*) portData.source();
  const int *fdims = field->lattice.dims();
  boost::multi_array_ref<short, 3>
    data((short *)field->lattice.dataPtr(),
	 boost::extents[fdims[2]][fdims[1]][fdims[0]]);

  // find number of projections
  int nprojections = 0;
  for (int i = 0; i < numberOfKeys; i++)
    if (imageKeyIds[i] == 0)
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
	   boost::extents[dims[2]][dims[1]][dims[0]]);
  // update angles
  int adims[3];
  adims[0] = nprojections;
  adims[1] = 1;
  adims[2] = 1;
  HxUniformScalarField3* new_angles =
    new HxUniformScalarField3(adims, McPrimType::mc_float);
  float *aptr = (float *)new_angles->lattice.dataPtr();
  int projectionIndex=0;
  int n_ibright=0;
  int n_fbright=0;
  int n_idark=0;
  int n_fdark=0;
  boost::multi_array<float, 2> i_bright(boost::extents[fdims[1]][fdims[0]]);
  boost::multi_array<float, 2> f_bright(boost::extents[fdims[1]][fdims[0]]);
  boost::multi_array<float, 2> i_dark(boost::extents[fdims[1]][fdims[0]]);
  boost::multi_array<float, 2> f_dark(boost::extents[fdims[1]][fdims[0]]);
  //initialize arrays
  for(int j=0;j<fdims[1];j++)
	  for(int i=0;i<fdims[0];i++)
		  i_bright[j][i]=f_bright[j][i]=i_dark[j][i]=f_dark[j][i]=0.0;
  for(int index=0;index<fdims[2];index++) // Loop through all the projection data
  {
	  //Check the Image key value
	  if(imageKeyIds[index]==0) //actual image
	  {
		  for(int j=0;j<fdims[1];j++)
			  for(int i=0;i<fdims[0];i++)
				pixels[projectionIndex][j][i] = data[index][j][i];
		  aptr[projectionIndex] = angles[index];
		  projectionIndex++;
	  }else if(imageKeyIds[index]==1){ //bright or open beam
		  if(index<fdims[2]/2){ //before experiment
			  for(int j=0;j<fdims[1];j++)
				  for(int i=0;i<fdims[0];i++)
					  i_bright[j][i] += data[index][j][i];
			  n_ibright++;
		  }else{ //after experiment
			  for(int j=0;j<fdims[1];j++)
				  for(int i=0;i<fdims[0];i++)
					  f_bright[j][i] += data[index][j][i];
			  n_fbright++;
		  }
	  }else if(imageKeyIds[index]==2){ //dark or closed beam
		  if(index<fdims[2]/2){ //before experiment
			  for(int j=0;j<fdims[1];j++)
				  for(int i=0;i<fdims[0];i++)
					  i_dark[j][i] += data[index][j][i];
			  n_idark++;
		  }else{ //after experiment
			  for(int j=0;j<fdims[1];j++)
				  for(int i=0;i<fdims[0];i++)
					  f_dark[j][i] += data[index][j][i];
			  n_fdark++;

		  }
	  }
  }

  //Average Brights field
  if(n_ibright>0){
	  for(int j=0;j<fdims[1];j++)
		  for(int i=0;i<fdims[0];i++)
			  i_bright[j][i] /= n_ibright;
  }
  if(n_fbright>0){
	  for(int j=0;j<fdims[1];j++)
		  for(int i=0;i<fdims[0];i++)
			  f_bright[j][i] /= n_fbright;
  }
  //Average Dark field
  if(n_idark>0){
	  for(int j=0;j<fdims[1];j++)
		  for(int i=0;i<fdims[0];i++)
			  i_dark[j][i] /= n_idark;
  }
  if(n_fdark>0){
	  for(int j=0;j<fdims[1];j++)
		  for(int i=0;i<fdims[0];i++)
			  f_dark[j][i] /= n_fdark;
  }
  // If no bright fields recorded then take the max intensity pixel
  if( n_ibright == 0 && n_fbright == 0)
  {
	  float max_pixel = pixels[0][0][0];
	  for (int k = 0; k < nprojections; k++) {
		  for (int j = 0; j < fdims[1]; j++) {
			  for (int i = 0; i < fdims[0]; i++) {
				  if (pixels[k][j][i] > max_pixel)
					  max_pixel = pixels[k][j][i];
			  }
		  }
	  }
	  //set the i_bright array and f_bright arrays for max pixel
	  for(int j=0;j<fdims[1];j++)
		  for(int i=0;i<fdims[0];i++)
			  i_bright[j][i]=f_bright[j][i]=max_pixel;
  }else if(n_ibright==0){
	  for(int j=0;j<fdims[1];j++)
		  for(int i=0;i<fdims[0];i++)
			  i_bright[j][i]=f_bright[j][i];
  }else if(n_fbright==0){
	  for(int j=0;j<fdims[1];j++)
		  for(int i=0;i<fdims[0];i++)
			  f_bright[j][i]=i_bright[j][i];
  }
  boost::multi_array<float, 2> dark(boost::extents[fdims[1]][fdims[0]]);
  boost::multi_array<float, 2> bright(boost::extents[fdims[1]][fdims[0]]);

  float minpixel = 0;
  float maxpixel = 0;
  //Normalise the data
  for(int index=0;index<projectionIndex;index++)
  {
	  // Based on fbp code, interpolate bright/dark
	  real w = (aptr[index] - aptr[0]) / (aptr[nprojections - 1] - aptr[0]);
	  for (int j = 0; j < fdims[1]; j++)
		  for (int i = 0; i < fdims[0]; i++)
			  dark[j][i] = i_dark[j][i] * (1.0 - w) + f_dark[j][i] * w;
	  for (int j = 0; j < fdims[1]; j++)
		  for (int i = 0; i < fdims[0]; i++)
			  bright[j][i] = i_bright[j][i] * (1.0 - w) + f_bright[j][i] * w;
	  // subtract dark from data/bright
	  // and clamp min data/bright value to 0.1
	  for (int j = 0; j < fdims[1]; j++) {
		  for (int i = 0; i < fdims[0]; i++) {
			  bright[j][i] -= dark[j][i];
			  if (bright[j][i] < 0.1)
				  bright[j][i] = 0.1;
		  }
	  }
	  for (int j = 0; j < fdims[1]; j++) {
		  for (int i = 0; i < fdims[0]; i++) {
			  pixels[index][j][i] -= dark[j][i];
			  if (pixels[index][j][i] < 0.1)
			  	  pixels[index][j][i] = 0.1;
		  }
	  }
	  // scale each data pixel by bright pixel
	  for (int j = 0; j < fdims[1]; j++)
		  for (int i = 0; i < fdims[0]; i++)
			  pixels[index][j][i] /= bright[j][i];

	  for (int j = 0; j < fdims[1]; j++) {
		  for (int i = 0; i < fdims[0]; i++) {
				 if(pixels[index][j][i] > maxpixel)
					 maxpixel = pixels[index][j][i];
				 if(pixels[index][j][i] < minpixel)
					 minpixel = pixels[index][j][i];
		  }
	  }	  
  }
  for(int index=0;index<projectionIndex;index++)
	  for(int j=0;j<fdims[1];j++)
		  for(int i=0;i<fdims[0];i++)
			  pixels[index][j][i] = (pixels[index][j][i]-minpixel)/(maxpixel-minpixel);
  // register result
  HxData::registerData(output, "projections");
  HxData::registerData(new_angles, "rotation angles");
  output->portMaster.connect(mcf);
  new_angles->portMaster.connect(mcf);
}
