/*
 *  Template of a compute module
 */

#ifndef CONE_BEAM_RECON_H
#define CONE_BEAM_RECON_H

#include <hxcore/HxCompModule.h>
#include <hxcore/HxPortDoIt.h>
#include <hxcore/HxPortIntSlider.h>
#include <hxcore/HxPortOnOff.h>
#include <hxcore/HxPortRadioBox.h>

#include "api.h"

class CCPIRECONSTRUCTION_API Cone_Beam_recon : public HxCompModule
{
    HX_HEADER(Cone_Beam_recon);

  public:

    Cone_Beam_recon();
    ~Cone_Beam_recon();

    HxPortDoIt portAction;
    HxConnection rotationAngle;
    // Todo - optional detector shift for some instruments
    //HxConnection horizontalShifts;
    //HxConnection verticalShifts;
    HxConnection pixelSize;
    HxConnection coneGeometry;
    HxPortRadioBox algorithm;
    HxPortIntSlider iterations;
    HxPortIntSlider resolution;
    HxPortOnOff beam_harden;

    virtual void compute();

 private:
  void run_reconstruction();
};

#endif // CONE_BEAM_RECON_H
