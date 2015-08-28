/**
 * This class is a implementation of CCPi User Application Interface in Avizo
 * Author: Mr. Srikanth Nagella
 * Date  : 15.05.2014
 */

#ifndef CCPIAVIZOUSERINTERFACE_H
#define CCPIAVIZOUSERINTERFACE_H

#include "CCPiUserApplicationInterface.h"

class CCPiAvizoUserInterface : public CCPiUserApplicationInterface
{
public:
	void LogMessage(std::string message);
	void SetStatusMessage(std::string message);
	void SetProgressValue(float value);
	bool isCancel();
};

#endif