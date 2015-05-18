/*
 *  Template of a compute module
 */

#include <QApplication>

#include <hxcore/HxMessage.h>
#include <hxfield/HxUniformScalarField3.h>

#include "NeXus_normalise.h"

HX_INIT_CLASS(NeXus_normalise,HxCompModule)

NeXus_normalise::NeXus_normalise() :
    HxCompModule(HxUniformScalarField3::getClassTypeId()),
    portAction(this,"action",QApplication::translate("NeXus_normalise", "Action"))
{
    portAction.setLabel(0,"DoIt");
}

NeXus_normalise::~NeXus_normalise()
{
}

void NeXus_normalise::compute()
{
    if (portAction.wasHit()) {
        theMsg->printf("do something\n");
    }
}
