/**
 * This class is a implementation of CCPi User Application Interface to output console
 * Author: Mr. Srikanth Nagella
 * Date  : 27.05.2014
 */

#ifndef CCPICONSOLEUSERINTERFACE_H
#define CCPICONSOLEUSERINTERFACE_H

#include "CCPiUserApplicationInterface.h"

class CCPiConsoleUserInterface : public CCPiUserApplicationInterface
{
public:
	void LogMessage(std::string message);
	void SetStatusMessage(std::string message);
	void SetProgressValue(float value);
	bool isCancel();
};

#endif