/*
 *  Template of a compute module
 */

#ifndef NEXUS_NORMALISE_H
#define NEXUS_NORMALISE_H

#include <hxcore/HxCompModule.h>
#include <hxcore/HxPortDoIt.h>

#include "api.h"
class HxMultiChannelField3;
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
  void normalise(HxMultiChannelField3 *);
  bool validateAndPopulateData(HxMultiChannelField3 *mcf);
  void normalise2();
  double *angles;
  int   *imageKeyIds;
  double *pixelSizeXY;
  int   numberOfKeys;
};

#endif // NEXUS_NORMALISE_H
