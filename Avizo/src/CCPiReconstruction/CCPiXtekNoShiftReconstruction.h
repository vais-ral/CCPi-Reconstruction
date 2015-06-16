/*
 *  Template of a compute module
 */

#ifndef CCPI_XTEK_NO_SHIFT_RECONSTRUCTION_H
#define CCPI_XTEK_NO_SHIFT_RECONSTRUCTION_H

#include <hxcore/HxCompModule.h>
#include <hxcore/HxPortDoIt.h>
#include <hxcore/HxPortIntSlider.h>
#include <hxcore/HxPortOnOff.h>
#include <hxcore/HxPortRadioBox.h>

#include "api.h"

class CCPIRECONSTRUCTION_API CCPiXtekNoShiftReconstruction : public HxCompModule
{
    HX_HEADER(CCPiXtekNoShiftReconstruction);

  public:

    CCPiXtekNoShiftReconstruction();
    ~CCPiXtekNoShiftReconstruction();

    HxPortDoIt portAction;
    HxPortRadioBox algorithm;
    HxPortIntSlider iterations;
    HxPortIntSlider resolution;
    HxPortOnOff beam_harden;

    virtual void compute();

 private:
  void run_reconstruction();
};

#endif // CCPI_XTEK_NO_SHIFT_RECONSTRUCTION_H
