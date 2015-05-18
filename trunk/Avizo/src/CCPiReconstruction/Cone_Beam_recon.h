/*
 *  Template of a compute module
 */

#ifndef CONE_BEAM_RECON_H
#define CONE_BEAM_RECON_H

#include <hxcore/HxCompModule.h>
#include <hxcore/HxPortDoIt.h>

#include "api.h"

class CCPIRECONSTRUCTION_API Cone_Beam_recon : public HxCompModule
{
    HX_HEADER(Cone_Beam_recon);

  public:

    Cone_Beam_recon();
    ~Cone_Beam_recon();

    HxPortDoIt portAction;

    virtual void compute();
};

#endif // CONE_BEAM_RECON_H
