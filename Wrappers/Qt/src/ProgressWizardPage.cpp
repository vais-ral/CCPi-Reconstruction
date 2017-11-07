#include "ProgressWizardPage.h"


ProgressWizardPage::ProgressWizardPage(void)
{
	bReconComplete = false;
}


ProgressWizardPage::~ProgressWizardPage(void)
{
}

bool ProgressWizardPage::isComplete() const
{
	return bReconComplete;
}

void  ProgressWizardPage::setReconComplete(bool value)
{
	bReconComplete = value;
	emit completeChanged();
}