/**
 * This class is a implementation of CCPi User Application Interface in Paraview
 * Author: Mr. Srikanth Nagella
 * Date  : 22.05.2014
 */

#ifndef CCPIPARAVIEWUSERINTERFACE_H
#define CCPIPARAVIEWUSERINTERFACE_H

#include "CCPiUserApplicationInterface.h"
#include "vtkAlgorithm.h"
#include "vtkSmartPointer.h"

class  CCPiParaviewUserInterface : public CCPiUserApplicationInterface
{
public:
    CCPiParaviewUserInterface(vtkAlgorithm *filter);
    void LogMessage(std::string message);
    void SetStatusMessage(std::string message);
    void SetProgressValue(float value);
    bool isCancel();
private:
    vtkSmartPointer<vtkAlgorithm> Filter;
};

#endif