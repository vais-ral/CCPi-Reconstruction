/*
 *  Template of a compute module
 */

#ifndef NEXUS_NORMALISE_H
#define NEXUS_NORMALISE_H

#include <hxcore/HxCompModule.h>
#include <hxcore/HxPortDoIt.h>

#include "api.h"

class CCPIRECONSTRUCTION_API NeXus_normalise : public HxCompModule
{
  HX_HEADER(NeXus_normalise);

 public:

  NeXus_normalise();
  ~NeXus_normalise();

  HxPortDoIt portAction;
  HxConnection rotationAngle;
  HxConnection imageKey;
  HxConnection pixelSize;

  virtual void compute();

 private:
  void normalise();
};

#endif // NEXUS_NORMALISE_H
