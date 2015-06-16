/**
 * This class is a filter to prepare the xtek data read from the files to prepare for reconstruction filter
 *
 * Author: Mr. Srikanth Nagella
 */
#ifndef CCPIXTEKAVIZOPREPAREFILTER_H
#define CCPIXTEKAVIZOPREPAREFILTER_H

#include <hxcore/HxCompModule.h>
#include <hxcore/HxPortIntSlider.h>
#include <hxcore/HxPortDoIt.h>
#include <hxfield/HxUniformScalarField3.h>
#include "api.h"

class CCPIRECONSTRUCTION_API CCPiXtekAvizoPrepareFilter : public HxCompModule
{
	HX_HEADER(CCPiXtekAvizoPrepareFilter);	

  public:
    CCPiXtekAvizoPrepareFilter();
    ~CCPiXtekAvizoPrepareFilter();

    /** Port providing a button to click to run the module */
    HxPortDoIt portAction;

    /** Perform the calculation in this module. Called by Avizo. */
    virtual void compute();

 private:
	int normalize(float *pixels, int numberOfProjections, int numberOfHorizontalPixels, int numberOfVerticalPixels, double whiteLevel,double scattering);
	int dataVerticalOffset;
	int numberOfVerticalPixels;
	int numberOfHorizontalPixels;
	float detectorPixelSize[2];
	float sourceToObject;
	float sourceToDetector;
	float maskRadius;
	double* angles;
	int numberOfProjections;
	HxUniformScalarField3* createOutput(HxUniformScalarField3 *field);
	void setParameters(HxUniformScalarField3* );
	void copyData(HxUniformScalarField3* src, HxUniformScalarField3* dst);

};

#endif
