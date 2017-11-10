#include "XradiaReader.h"
#include "api.h"
#include <hxcore/HxData.h>
#include <hxfield/HxUniformScalarField3.h>
#include <amiramesh/HxParamBundle.h>
#include "CCPiAvizoUserInterface.h"
#if AVIZO_UNSUPPORTED_INTERNAL
#include <hxcore/internal/HxWorkArea.h>
#include <hxcore/internal/HxThread.h>
#else
#include <hxcore/HxWorkArea.h>
#include <hxcore/HxThread.h>
#endif
#include <QApplication.h>
#include <boost/filesystem.hpp>

void setParameters(HxUniformScalarField3 *field, CCPi::XradiaReader reader);

CCPI_API
int CCPiXradiaAvizoReader(const char* filename)
{
    theWorkArea->startWorking(
     QApplication::translate("CCPiXradiaAvizoReader", "Reading xradia data"));   
    CCPiAvizoUserInterface *msg = new CCPiAvizoUserInterface();
	CCPi::XradiaReader reader(filename,msg);
	int newDims[3];
	newDims[0] = reader.getImageWidth();
	newDims[1] = reader.getImageHeight();
	newDims[2] = reader.getNumberOfImages();
#if AVIZO_UNSUPPORTED_INTERNAL
	HxUniformScalarField3 *field = new HxUniformScalarField3(newDims, McPrimType::MC_UINT16);
	uint16_t *data = (uint16_t*)reader.getImageDataInShort()->data();
	uint16_t *rawout=(uint16_t *)field->lattice().dataPtr();
	unsigned long long totalsize = 1;
	for(int i=0;i<3;i++)
		totalsize *= newDims[i];

	for(unsigned long long idx=0;idx<totalsize;idx++)
		rawout[idx]=data[idx];

	HxUniformCoord3 * coords = (HxUniformCoord3 *)field->lattice().coords();

	McBox3f bx;
	bx = coords->getBoundingBox();
#else
	HxUniformScalarField3 *field = new HxUniformScalarField3(newDims, McPrimType::mc_uint16);
	uint16_t *data = (uint16_t*)reader.getImageDataInShort()->data();
	uint16_t *rawout = (uint16_t *)field->lattice.dataPtr();
	unsigned long long totalsize = 1;
	for (int i = 0; i<3; i++)
		totalsize *= newDims[i];

	for (unsigned long long idx = 0; idx<totalsize; idx++)
		rawout[idx] = data[idx];

	HxUniformCoord3 * coords = (HxUniformCoord3 *)field->lattice.coords();

	float * bx;
	bx=coords->bbox();
#endif
	bx[0]=0;
	if(newDims[0]!=1)
		bx[1]=(float)newDims[0]-1;
	else
		bx[1]=1;
	bx[2]=0;
	if(newDims[1]!=1)
		bx[3]=(float)newDims[1]-1;
	else
		bx[3]=1;
	bx[4]=0;
	if(newDims[2]!=1)
		bx[5]=(float)newDims[2]-1;
	else
		bx[5] = 1;
#if AVIZO_UNSUPPORTED_INTERNAL
	coords->setBoundingBox(bx);
#endif
	HxData::registerData(field, "ImageData");
	setParameters(field, reader);
	theWorkArea->stopWorking();
	return 1;
}

void setParameters(HxUniformScalarField3 *field, CCPi::XradiaReader reader)
{
	//field->parameters.set("VoxelSize", 3, reader.getVoxelSize());
	//field->parameters.set("Offset", 3, reader.getOffset());
	field->parameters.set("SourceToObject", -1*reader.getSourceToObject());
	field->parameters.set("SourceToDetector", reader.getDetectorToObject()-reader.getSourceToObject());
	//Copy angles
	double *angles = new double[reader.getNumberOfImages()];
	std::vector<float> rangles = reader.getAngles();
	for(int i=0;i<reader.getNumberOfImages();i++)
		angles[i] = rangles[i];
	field->parameters.set("Angles", reader.getNumberOfImages(), angles);
	//field->parameters.set("DetectorPixelSize", 2, reader.getDetectorPixelSize());
	////field->parameters.set("InitialAngle", reader.getInitialAngle());
	//field->parameters.set("WhiteLevel", reader.getWhiteLevel());
	//field->parameters.set("Scattering", reader.getScattering());
	//field->parameters.set("CoefX", 5, reader.getCoefX());
	//field->parameters.set("Scale", reader.getScale());
	//field->parameters.set("MaskRadius", reader.getMaskRadius());
}
