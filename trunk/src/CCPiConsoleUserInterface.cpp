#include "CCPiConsoleUserInterface.h"
#include <iostream>

void CCPiConsoleUserInterface::LogMessage(std::string message)
{
	std::cout << message.c_str() <<std::endl;
}

void CCPiConsoleUserInterface::SetStatusMessage(std::string message)
{
	std::cout<< "Status: "<<message.c_str()<<std::endl;
}

void CCPiConsoleUserInterface::SetProgressValue(float value)
{
	std::cout<< "Progress: "<<value<<std::endl;
}

bool CCPiConsoleUserInterface::isCancel()
{
	return false;
}
