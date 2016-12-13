#ifndef CCPI_SAVINGWIZARD_PAGE_H
#define CCPI_SAVINGWIZARD_PAGE_H

#include <QWizardPage>
class SavingWizardPage :
	public QWizardPage
{
public:
	SavingWizardPage(void);
	~SavingWizardPage(void);
	bool isComplete() const;
	void setSavingComplete(bool value);
private:
	bool bSavingComplete;
};

#endif
