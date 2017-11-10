#include "CCPiAvizoUserInterface.h"

#include <hxcore/HxMessage.h>
#include <hxcore/internal/HxWorkArea.h>


void CCPiAvizoUserInterface::LogMessage(std::string message)
{
	theMsg->stream() << message.c_str() <<std::endl;
}

void CCPiAvizoUserInterface::SetStatusMessage(std::string message)
{
	#if QT_VERSION >= 0x050000
		theWorkArea->setProgressInfo(QString::fromStdString(message));
	#else
		theWorkArea->setProgressInfo(QString(message.c_str()));
	#endif
}

void CCPiAvizoUserInterface::SetProgressValue(float value)
{
	theWorkArea->setProgressValue(value);
}

bool CCPiAvizoUserInterface::isCancel()
{
	return theWorkArea->wasInterrupted();
}
