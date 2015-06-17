#include <QApplication>
#include <hxcore/HxMessage.h>
#include <hxcore/HxWorkArea.h>
#include <hxfield/HxUniformScalarField3.h>
#include "CCPiXtekAvizoPrepareFilter.h"

HX_INIT_CLASS(CCPiXtekAvizoPrepareFilter,HxCompModule)

CCPiXtekAvizoPrepareFilter::CCPiXtekAvizoPrepareFilter() :
    HxCompModule(HxUniformScalarField3::getClassTypeId()),
    portAction(this,"action",QApplication::translate("CCPiXtekAvizoPrepareFilter", "Action"))
{
    portAction.setLabel(0,"DoIt");
}

CCPiXtekAvizoPrepareFilter::~CCPiXtekAvizoPrepareFilter()
{
}

void CCPiXtekAvizoPrepareFilter::compute()
{
	// Check whether the action port button was clicked
    if (!portAction.wasHit()) return;

	// Access the input data object. The member portData, which is of type
    // HxConnection, is inherited from HxModule.
    HxUniformScalarField3 *field = (HxUniformScalarField3*) portData.source();
    // Check input data
	if (field->primType() != McPrimType::mc_float || field->parameters.find("VoxelSize",3) == NULL) {
        theMsg->stream() << "This module only works on a Xtek data" <<std::endl;
        return;
    }
	numberOfVerticalPixels = field->lattice.dims()[1];
	numberOfHorizontalPixels = field->lattice.dims()[0];
	numberOfProjections = field->lattice.dims()[2];
	angles = new double[numberOfProjections];
	double pixelSize[2];
	field->parameters.findReal("DetectorPixelSize",2,pixelSize);
	detectorPixelSize[0] = pixelSize[0];
	detectorPixelSize[1] = pixelSize[1];

	//source to Detector and object
	field->parameters.findReal("SourceToDetector", sourceToDetector);
	field->parameters.findReal("SourceToObject", sourceToObject);
	field->parameters.findReal("MaskRadius", maskRadius);

    // Create an output with same size as input. Data type will be unsigned char
    // as we produce a labelled image.
	HxUniformScalarField3 *output = createOutput(field);
    // Output shall have same bounding box as input
    output->coords()->setBoundingBox(field->bbox());
	//Set the output parameters
	setParameters(output);
	//Copy the data
	copyData(field, output);
	//normalize data
	double whiteLevel;
	double scattering;
	field->parameters.findReal("WhiteLevel",whiteLevel);
	field->parameters.findReal("Scattering",scattering);
	normalize((float*)output->lattice.dataPtr(), numberOfProjections, numberOfHorizontalPixels,numberOfVerticalPixels, whiteLevel, scattering);
	// Add to the workspace
	setResult(output);
    // Stop progress bar
    theWorkArea->stopWorking();
}


/**
 * Create an output with same size as input. The output image is essentially
 * modified version of the input data.
 * Check if there is a result that can be re-used. If so check type and size
 * match current input.
 * @param field The input field
 * @return The output to use for showing data in application
 */
HxUniformScalarField3* CCPiXtekAvizoPrepareFilter::createOutput(HxUniformScalarField3 *field)
{
    // Check if there is a result which we can reuse
    HxUniformScalarField3 *output = (HxUniformScalarField3*) getResult();
    
    // Check for proper type
    if ( output && !output->isOfType(HxUniformScalarField3::getClassTypeId()) )
        output = NULL;
    
    // Check if size and primitive type still match current input
    const int *dims = field->lattice.dims();
    if (output) {
        const int *outdims = output->lattice.dims();
        if ( dims[0] != outdims[0] || dims[1] != outdims[1] ||
            dims[2] != outdims[2] || output->primType() != McPrimType::mc_float )
            
            output = NULL;
    }
    
    // If necessary create a new result data set
    if (!output) {
        output = new HxUniformScalarField3(dims, McPrimType::mc_float);
        output->composeLabel(field->getName(), "result");
    }

    return output;
}

/**
 * Copies the data from src to dst.
 * @param src source data
 * @param dst destination data
 */
void CCPiXtekAvizoPrepareFilter::copyData(HxUniformScalarField3* src, HxUniformScalarField3* dst)
{
	float *dstPixels = (float *)dst->lattice.dataPtr();
	float *srcPixels = (float *)src->lattice.dataPtr();
	//Change the data image order
	long long angval = 1;
	long long vval = 1;
	for(int ang=0;ang<numberOfProjections;ang++){
		angval = ang * numberOfHorizontalPixels * numberOfVerticalPixels;
		for (int v = 0; v < numberOfVerticalPixels; v++) {
			vval = v*numberOfHorizontalPixels;
			for (int h = 0; h < numberOfHorizontalPixels; h++) {
				dstPixels[angval+vval+h] = srcPixels[angval+(numberOfVerticalPixels - v - 1) * numberOfHorizontalPixels + h];
			}
		}
	}
}

/**
 * Sets/copies the default parameter values for the field
 * @param the input scalar field for which the parameter are set
 */
void CCPiXtekAvizoPrepareFilter::setParameters(HxUniformScalarField3* field)
{
	field->parameters.set("Angles", numberOfProjections, angles);
	field->parameters.set("DetectorPixelSize", 2, detectorPixelSize);
	field->parameters.set("SourceToObject", sourceToObject);
	field->parameters.set("SourceToDetector", sourceToDetector);
	field->parameters.set("MaskRadius", maskRadius);
}

/**
 * This method normalises the pixel data using white level and scattering. 
 * @param pixels pixel data
 * @param iNumberOfProjections number of projection 
 * @param iNumberOfHorizontalPixels number of horizontal pixels
 * @param iNumberOfVerticalPixels number of vertical pixels
 * @param whiteLevel white level value for the experiment
 * @param scattering scattering setting for the experiment.
 * @return whether the normalisation of pixel data was success full.
 */
int CCPiXtekAvizoPrepareFilter::normalize(float *pixels, int iNumberOfProjections, int iNumberOfHorizontalPixels, int iNumberOfVerticalPixels, double whiteLevel,double scattering)
{
	float max_v = 0.0;
	bool fail = false;
	unsigned long long index = 0;
	//Check that the pixel value is always less than whiteLevel value
	for (int i = 0; i < iNumberOfProjections; i++) {
		for (int k = 0; k < iNumberOfVerticalPixels; k++) {
			for (int j = 0; j < iNumberOfHorizontalPixels; j++) {
				if (max_v < pixels[index])
					max_v = pixels[index];
				if (max_v > whiteLevel)
					fail = true;
				index++;
			}
		}
	}
	if (fail) {
		theMsg->stream() << "Values exceed white level"<<std::endl;
//		return false;
	} else
	  max_v = whiteLevel;
	max_v -= whiteLevel * scattering / float(100.0); 
	for (index = 0; index < iNumberOfProjections*iNumberOfHorizontalPixels*(unsigned long long)iNumberOfVerticalPixels; index++) {
		pixels[index] = (pixels[index]
			      - whiteLevel * scattering / float(100.0)) / max_v;
	}
	return true;
}
