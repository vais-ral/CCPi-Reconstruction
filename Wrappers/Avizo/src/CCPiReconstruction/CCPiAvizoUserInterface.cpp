#include "CCPiAvizoUserInterface.h"

#include <hxcore/HxMessage.h>
#include <hxcore/internal/HxWorkArea.h>


void CCPiAvizoUserInterface::LogMessage(std::string message)
{
	theMsg->stream() << message.c_str() <<std::endl;
}

void CCPiAvizoUserInterface::SetStatusMessage(std::string message)
{
	theWorkArea->setProgressInfo(QString(message.c_str()));
}

void CCPiAvizoUserInterface::SetProgressValue(float value)
{
	theWorkArea->setProgressValue(value);
}

bool CCPiAvizoUserInterface::isCancel()
{
	return theWorkArea->wasInterrupted();
}
