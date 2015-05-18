/*
 *  Template of a compute module
 */

#include <QApplication>

#include <hxcore/HxMessage.h>
#include <hxfield/HxUniformScalarField3.h>

#include "Cone_Beam_recon.h"

HX_INIT_CLASS(Cone_Beam_recon,HxCompModule)

Cone_Beam_recon::Cone_Beam_recon() :
    HxCompModule(HxUniformScalarField3::getClassTypeId()),
    portAction(this,"action",QApplication::translate("Cone_Beam_recon", "Action"))
{
    portAction.setLabel(0,"DoIt");
}

Cone_Beam_recon::~Cone_Beam_recon()
{
}

void Cone_Beam_recon::compute()
{
    if (portAction.wasHit()) {
        theMsg->printf("do something\n");
    }
}
