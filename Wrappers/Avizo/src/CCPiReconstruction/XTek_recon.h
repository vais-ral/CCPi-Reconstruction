/*
 *  Template of a compute module
 */

#ifndef XTEK_RECON_H
#define XTEK_RECON_H

#include <hxcore/HxCompModule.h>
#include <hxcore/HxPortDoIt.h>
#include <hxcore/HxPortIntSlider.h>
#include <hxcore/HxPortOnOff.h>
#include <hxcore/HxPortRadioBox.h>
#include <hxcore/HxFileDialog.h>

#include "api.h"

class CCPI_API XTek_recon : public HxCompModule
{
    HX_HEADER(XTek_recon);

  public:

    HxPortDoIt portAction;
    HxPortRadioBox algorithm;
    HxPortIntSlider iterations;
    HxPortIntSlider resolution;
    HxPortOnOff beam_harden;
    HxFileDialog browser;

    virtual void compute();

 private:
    void run_reconstruction(const std::string filename);
};

#endif // XTEK_RECON_H
