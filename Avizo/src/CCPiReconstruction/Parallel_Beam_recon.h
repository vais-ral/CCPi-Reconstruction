/*
 *  Template of a compute module
 */

#ifndef PARALLEL_BEAM_RECON_H
#define PARALLEL_BEAM_RECON_H

#include <hxcore/HxCompModule.h>
#include <hxcore/HxPortDoIt.h>
#include <hxcore/HxPortIntSlider.h>
#include <hxcore/HxPortOnOff.h>
#include <hxcore/HxPortRadioBox.h>
#include <hxcore/HxPortFloatTextN.h>

#include "api.h"

class CCPIRECONSTRUCTION_API Parallel_Beam_recon : public HxCompModule
{
  HX_HEADER(Parallel_Beam_recon);

 public:

  Parallel_Beam_recon();
  ~Parallel_Beam_recon();

  HxPortDoIt portAction;
  HxConnection rotationAngle;
  HxConnection pixelSize;
  HxPortRadioBox algorithm;
  HxPortIntSlider iterations;
  HxPortIntSlider resolution;
  HxPortFloatTextN rotationCentre;
  HxPortOnOff beam_harden;

  virtual void compute();

 private:
  void run_reconstruction();
};

#endif // PARALLEL_BEAM_RECON_H
