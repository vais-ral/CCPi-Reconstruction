#include "XtekReader.h"
#include "api.h"
#include <hxcore/HxData.h>
#include <hxfield/HxUniformScalarField3.h>
#include <amiramesh/HxParamBundle.h>
#include "CCPiAvizoUserInterface.h"
#include <hxcore/internal/HxWorkArea.h>
#include <QApplication>
#include <hxcore/internal/HxThread.h>
#include <boost/filesystem.hpp>

void setParameters(HxUniformScalarField3 *field, XtekReader reader);

CCPI_API
int CCPiXTekAvizoReader(const char* filename)
{
//    theWorkArea->startWorking(
//     QApplication::translate("CCPiXTekAvizoReader", "Reading xtekct data"));   
	 CCPiAvizoUserInterface *msg = new CCPiAvizoUserInterface();
	XtekReader reader(filename,msg);
	int newDims[3];
	newDims[0] = reader.getImageWidth();
	newDims[1] = reader.getImageHeight();
	newDims[2] = reader.getNumberOfProjections();

#if AVIZO_UNSUPPORTED_INTERNAL
	HxUniformScalarField3 *field = new HxUniformScalarField3(newDims, McPrimType::MC_UINT16);
	uint16_t *data = (uint16_t*)reader.getImageData();
	uint16_t *rawout = (uint16_t *)field->lattice().dataPtr();
	unsigned long long totalsize = 1;
	for (int i = 0; i<3; i++)
		totalsize *= newDims[i];

	for (unsigned long long idx = 0; idx<totalsize; idx++)
		rawout[idx] = data[idx];

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
	bx = coords->bbox();
#endif
	bx[0] = 0;
	if (newDims[0] != 1)
		bx[1] = (float)newDims[0] - 1;
	else
		bx[1] = 1;
	bx[2] = 0;
	if (newDims[1] != 1)
		bx[3] = (float)newDims[1] - 1;
	else
		bx[3] = 1;
	bx[4] = 0;
	if (newDims[2] != 1)
		bx[5] = (float)newDims[2] - 1;
	else
		bx[5] = 1;
#if AVIZO_UNSUPPORTED_INTERNAL
	coords->setBoundingBox(bx);
#endif
#if QT_VERSION >= 0x050000
	HxData::registerData(field, QString::fromStdString(reader.getName()));
#else
	HxData::registerData(field, reader.getName().c_str());
#endif
	setParameters(field, reader);
	theWorkArea->stopWorking();
	return 1;
}

void setParameters(HxUniformScalarField3 *field, XtekReader reader)
{
	field->parameters.set("VoxelSize", 3, reader.getVoxelSize());
	field->parameters.set("Offset", 3, reader.getOffset());
	field->parameters.set("SourceToObject", reader.getSourceToObject());
	field->parameters.set("SourceToDetector", reader.getSourceToDetector());
	field->parameters.set("Angles", reader.getNumberOfProjections(), reader.getAngles());
	field->parameters.set("DetectorPixelSize", 2, reader.getDetectorPixelSize());
	//field->parameters.set("InitialAngle", reader.getInitialAngle());
	field->parameters.set("WhiteLevel", reader.getWhiteLevel());
	field->parameters.set("Scattering", reader.getScattering());
	field->parameters.set("CoefX", 5, reader.getCoefX());
	field->parameters.set("Scale", reader.getScale());
	field->parameters.set("MaskRadius", reader.getMaskRadius());
	
}
