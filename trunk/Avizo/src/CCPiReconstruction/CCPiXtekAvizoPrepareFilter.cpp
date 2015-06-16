#include <QApplication>
#include <hxcore/HxMessage.h>
#include <hxcore/HxWorkArea.h>
#include <hxfield/HxUniformScalarField3.h>
#include "CCPiXtekAvizoPrepareFilter.h"
#include "base_types.hpp"

HX_INIT_CLASS(CCPiXtekAvizoPrepareFilter,HxCompModule)

CCPiXtekAvizoPrepareFilter::CCPiXtekAvizoPrepareFilter() :
    HxCompModule(HxUniformScalarField3::getClassTypeId()),
    portAction(this,"action",QApplication::translate("CCPiXtekAvizoPrepareFilter", "Action")),
    resolution(this, "resolution", QApplication::translate("XTek_recon", "Pixels per Voxel"))
{
    portAction.setLabel(0,"DoIt");
	resolution.setMinMax(1, 8);
	resolution.setValue(1);
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

	numberOfVerticalPixels = field->lattice.dims()[0];
	numberOfHorizontalPixels = field->lattice.dims()[1];
	numberOfProjections = field->lattice.dims()[2];
	///// original xtek.cpp translation
	int npix = calculateVerticalAlignment(numberOfVerticalPixels, resolution.getValue(), true);
	int diff = npix - numberOfVerticalPixels;
	dataVerticalOffset = diff / 2;
	// try and check that midpoint is in aligned middle of pixels?
	// if n is an even number then they must be, mustn't they?
	numberOfVerticalPixels = npix;
	verticalPixels = new double[numberOfVerticalPixels];
	horizontalPixels = new double[numberOfHorizontalPixels];
	double pixelSize[2];
	field->parameters.findReal("DetectorPixelSize",2,pixelSize);
	detectorPixelSize[0] = pixelSize[0];
	detectorPixelSize[1] = pixelSize[1];
	real pixel_base = -((numberOfHorizontalPixels - 1) * pixelSize[0] / real(2.0));
	for (int i = 0; i < numberOfHorizontalPixels; i++)
		horizontalPixels[i] = pixel_base + real(i) * pixelSize[0];
	pixel_base = -((numberOfVerticalPixels - 1) * pixelSize[1] / real(2.0));
	for (int i = 0; i < numberOfVerticalPixels; i++)
		horizontalPixels[i] = pixel_base + real(i) * pixelSize[1];

	//convert angles to radians
	angles = new double[numberOfProjections];
	field->parameters.findReal("Angles", numberOfProjections, angles);
	for(int i=0;i<numberOfProjections;i++)
		angles[i] = angles[i] *  real(M_PI)/real(180.0);

	//source to Detector and object
	field->parameters.findReal("SourceToDetector", sourceToDetector);
	field->parameters.findReal("SourceToObject", sourceToObject);

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
	long long hval = 1;
	for(int ang=0;ang<numberOfProjections;ang++){
		angval = ang * numberOfHorizontalPixels * numberOfVerticalPixels;
		for (int h = 0; h < numberOfHorizontalPixels; h++) {
			hval = h*numberOfVerticalPixels;
			for (int v = 0; v < numberOfVerticalPixels; v++) {
				dstPixels[angval+hval+v] = srcPixels[angval+(numberOfVerticalPixels - v - 1) * numberOfHorizontalPixels + h];
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
	field->parameters.set("Resolution", resolution.getValueToInt());
}

/**
 * Method calculates the vertical alignment
 * @param n number of pixels
 * @param pix_per_vox number of pixels per voxel
 * @param cone whether the beam type is cone
 * @return the number of pixels after vertical alignment
 */
int CCPiXtekAvizoPrepareFilter::calculateVerticalAlignment(const int n, const int pix_per_vox, const bool cone)
{
	int nvox = n / pix_per_vox;
	if (n % pix_per_vox != 0)
		nvox++;
	int align = aligned_allocator<voxel_type>::alignment / sizeof(voxel_type);
	int vox_blocks = nvox / align;
	if (nvox % align != 0)
		vox_blocks++;
	// make sure the cone can split into 2 equal aligned halves
	if (cone) {
		if (vox_blocks %2 != 0)
			vox_blocks++;
	}
	int npix = vox_blocks * align * pix_per_vox;
	return npix;
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
	real max_v = 0.0;
	bool fail = false;
	unsigned long long index = 0;
	unsigned long long tmpIndexI = 0;
	unsigned long long tmpIndexJ = 0;
	//Check that the pixel value is always less than whiteLevel value
	for (int i = 0; i < iNumberOfProjections; i++) {
		tmpIndexI = i * iNumberOfHorizontalPixels * iNumberOfVerticalPixels;
		for (int j = 0; j < iNumberOfHorizontalPixels; j++) {
			tmpIndexJ = j*iNumberOfHorizontalPixels;
			for (int k = 0; k < iNumberOfVerticalPixels; k++) {
				index = tmpIndexI + tmpIndexJ + k;
				if (max_v < pixels[index])
					max_v = pixels[index];
				if (max_v > whiteLevel)
					fail = true;
			}
		}
	}
	if (fail) {
		theMsg->stream() << "Values exceed white level"<<std::endl;
//		return false;
	}
	max_v = whiteLevel; 
	// scale and take -ve log, due to exponential extinction in sample.
	for (index = 0; index < iNumberOfProjections*iNumberOfHorizontalPixels*(unsigned long long)iNumberOfVerticalPixels; index++) {
				pixels[index] -= whiteLevel * scattering / real(100.0);
				if (pixels[index] < real(1.0))
					pixels[index] = - std::log(real(0.00001) / max_v);
				else
					pixels[index] = - std::log(pixels[index] / max_v);
	}
	return true;
}