/**
 * This class is a filter to prepare the xradia data read from the files to prepare for reconstruction filter
 *
 * Author: Mr. Srikanth Nagella
 */
#ifndef CCPIXRADIAAVIZOPREPAREFILTER_H
#define CCPIXRADIAAVIZOPREPAREFILTER_H

#include <hxcore/HxCompModule.h>
#include <hxcore/HxPortIntSlider.h>
#include <hxcore/HxPortDoIt.h>
#include <hxfield/HxUniformScalarField3.h>
#include "api.h"

class CCPIRECONSTRUCTION_API CCPiXradiaAvizoPrepareFilter : public HxCompModule
{
	HX_HEADER(CCPiXradiaAvizoPrepareFilter);	

  public:

    /** Port providing a button to click to run the module */
    HxPortDoIt portAction;

    /** Perform the calculation in this module. Called by Avizo. */
    virtual void compute();

 private:
	int dataVerticalOffset;
	int numberOfVerticalPixels;
	int numberOfHorizontalPixels;
	float detectorPixelSize[2];
	float sourceToObject;
	float sourceToDetector;
	float* angles;
	int numberOfProjections;
	HxUniformScalarField3* createOutput(HxUniformScalarField3 *field);
	void setParameters(HxUniformScalarField3* );
	void copyData(HxUniformScalarField3* src, HxUniformScalarField3* dst);

};

#endif
