#include <QApplication>
#include <hxcore/HxMessage.h>
#include <hxcore/internal/HxWorkArea.h>
#include <hxfield/HxUniformScalarField3.h>
#include "CCPiXradiaAvizoPrepareFilter.h"

HX_INIT_CLASS(CCPiXradiaAvizoPrepareFilter,HxCompModule)

CCPiXradiaAvizoPrepareFilter::CCPiXradiaAvizoPrepareFilter() :
    HxCompModule(HxUniformScalarField3::getClassTypeId()),
    portAction(this,"action",QApplication::translate("CCPiXradiaAvizoPrepareFilter", "Action"))
{
    portAction.setLabel(0,"DoIt");
}

CCPiXradiaAvizoPrepareFilter::~CCPiXradiaAvizoPrepareFilter()
{
}

void CCPiXradiaAvizoPrepareFilter::compute()
{
	// Check whether the action port button was clicked
    if (!portAction.wasHit()) return;

	// Access the input data object. The member portData, which is of type
    // HxConnection, is inherited from HxModule.
    HxUniformScalarField3 *field = (HxUniformScalarField3*) portData.getSource();
    // Check input data
	if (field->primType() != McPrimType::MC_UINT16 || field->parameters.find("SourceToDetector",1) == NULL) {
        theMsg->stream() << "This module only works on a Xradia data" <<std::endl;
        return;
    }
    theWorkArea->startWorking(
     QApplication::translate("CCPiXradiaAvizoPrepareFilter", "Preparing data"));   
	numberOfVerticalPixels = field->lattice().getDims()[1];
	numberOfHorizontalPixels = field->lattice().getDims()[0];
	numberOfProjections = field->lattice().getDims()[2];
	angles = new float[numberOfProjections];
	field->parameters.findReal("Angles", numberOfProjections, angles);
	  
	double pixelSize[2];
//	field->parameters.findReal("DetectorPixelSize",2,pixelSize);
	detectorPixelSize[0] = 1.0;
	detectorPixelSize[1] = 1.0;

	//source to Detector and object
	field->parameters.findReal("SourceToDetector", sourceToDetector);
	field->parameters.findReal("SourceToObject", sourceToObject);
	theWorkArea->setProgressInfo("Allocating memory");
	theWorkArea->setProgressValue(0.1);
    // Create an output with same size as input. Data type will be unsigned char
    // as we produce a labelled image.
	HxUniformScalarField3 *output = createOutput(field);
	theWorkArea->setProgressInfo("Copying data");
	theWorkArea->setProgressValue(0.4);
    // Output shall have same bounding box as input
    output->coords()->setBoundingBox(field->getBoundingBox());
	//Set the output parameters
	setParameters(output);
	//Copy the data
	copyData(field, output);
	//TODO:: Normalize the data data. don't know how it's done in Xradia
	theWorkArea->setProgressValue(1.0);
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
HxUniformScalarField3* CCPiXradiaAvizoPrepareFilter::createOutput(HxUniformScalarField3 *field)
{
    // Check if there is a result which we can reuse
    HxUniformScalarField3 *output = (HxUniformScalarField3*) getResult();
    
    // Check for proper type
    if ( output && !output->isOfType(HxUniformScalarField3::getClassTypeId()) )
        output = NULL;
    
    // Check if size and primitive type still match current input
    McDim3l dims = field->lattice().getDims();
    if (output) {
        McDim3l outdims = output->lattice().getDims();
        if ( dims[0] != outdims[0] || dims[1] != outdims[1] ||
			dims[2] != outdims[2] || output->primType() != McPrimType::MC_FLOAT )
            
            output = NULL;
    }
    
    // If necessary create a new result data set
    if (!output) {
        output = new HxUniformScalarField3(dims, McPrimType::MC_FLOAT);
        output->composeLabel(field->getName(), "result");
    }

    return output;
}

/**
 * Copies the data from src to dst.
 * @param src source data
 * @param dst destination data
 */
void CCPiXradiaAvizoPrepareFilter::copyData(HxUniformScalarField3* src, HxUniformScalarField3* dst)
{
	float *dstPixels = (float *)dst->lattice().dataPtr();
	uint16_t *srcPixels = (uint16_t *)src->lattice().dataPtr();
	//Change the data image order
	unsigned long long angval = 1;
	unsigned long long vval = 1;
	for(unsigned long long ang=0;ang<numberOfProjections;ang++){
		angval = ang * numberOfHorizontalPixels * numberOfVerticalPixels;
		for (unsigned long long v = 0; v < numberOfVerticalPixels; v++) {
			vval = v*numberOfHorizontalPixels;
			for (unsigned long long h = 0; h < numberOfHorizontalPixels; h++) {
				dstPixels[angval+vval+h] = srcPixels[angval+(numberOfVerticalPixels - v - 1) * numberOfHorizontalPixels + h];
			}
		}
	}
}

/**
 * Sets/copies the default parameter values for the field
 * @param the input scalar field for which the parameter are set
 */
void CCPiXradiaAvizoPrepareFilter::setParameters(HxUniformScalarField3* field)
{
	field->parameters.set("Angles", numberOfProjections, angles);
	field->parameters.set("DetectorPixelSize", 2, detectorPixelSize);
	field->parameters.set("SourceToObject", sourceToObject);
	field->parameters.set("SourceToDetector", sourceToDetector);
	field->parameters.set("MaskRadius", 1.0);//TODO:: Correct this as there is no mask radius in Xradia
}

