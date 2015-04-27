/*
 *  Template of a compute module
 */

#ifndef CGLS_RECON_H
#define CGLS_RECON_H

#include <hxcore/HxCompModule.h>
#include <hxcore/HxPortDoIt.h>
#include <hxcore/HxPortIntSlider.h>
#include <hxcore/HxPortOnOff.h>

#include "api.h"

class CCPI_API CGLS_recon : public HxCompModule
{
  HX_HEADER(CGLS_recon);

 public:

  CGLS_recon();
  ~CGLS_recon();
  /** Connection to masking data */
  HxConnection rotationAngle;
  HxConnection pixelSize;
  HxConnection imageKey;
  HxPortDoIt portAction;
  HxPortIntSlider iterations;
  HxPortIntSlider resolution;
  HxPortOnOff beam_harden;

  virtual void compute();

 private:
  void run_cgls();
};

#endif // CGLS_RECON_H
