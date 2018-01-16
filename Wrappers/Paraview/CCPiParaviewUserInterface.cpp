#include "CCPiParaviewUserInterface.h"
#include <string>

CCPiParaviewUserInterface::CCPiParaviewUserInterface(vtkAlgorithm *filter)
{
    Filter = filter;
}

void CCPiParaviewUserInterface::LogMessage(std::string message)
{
    vtkGenericWarningMacro(<<message.c_str());
}

void CCPiParaviewUserInterface::SetStatusMessage(std::string message)
{
    Filter->SetProgressText(message.c_str());
}

void CCPiParaviewUserInterface::SetProgressValue(float value)
{
    Filter->UpdateProgress(value);
}

bool CCPiParaviewUserInterface::isCancel()
{
    return false;
}