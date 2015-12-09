#include "SavingWizardPage.h"


SavingWizardPage::SavingWizardPage(void)
{
	bSavingComplete = false;
}


SavingWizardPage::~SavingWizardPage(void)
{
}

bool SavingWizardPage::isComplete() const
{
	return bSavingComplete;
}

void  SavingWizardPage::setSavingComplete(bool value)
{
	bSavingComplete = value;
	emit completeChanged();
}