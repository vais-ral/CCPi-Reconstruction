/**
 * This is an interface between CCPi core algorithms and External application specific calls.
 * such as Avizo, etc.
 * Author: Mr. Srikanth Nagella
 * Date  : 14.05.2014
 */

#ifndef CCPIUSERAPPLICATIONINTERFACE_H
#define CCPIUSERAPPLICATIONINTERFACE_H

#include <string>

class CCPiUserApplicationInterface
{
public:
	/* 
	 * Message to logging unit such as console window or files
	 * @message : string to be send to logger
	 */
	virtual void LogMessage(std::string message)=0;
	/* Current status update in the user interface
	 * @message : current status string to be set in the user interface
	 */
	virtual void SetStatusMessage(std::string message)=0;
	/* Progress of the current running algorithm 
	 * @value : set the progress between 0.0(started) - 1.0 (completed)
	 */
	virtual void SetProgressValue(float value)=0;
	/* Whether the cancel/stop button is pressed in the user interface so that
	   algorithm can exit immediately without completing*/
	virtual bool isCancel()=0;
};

#endif
