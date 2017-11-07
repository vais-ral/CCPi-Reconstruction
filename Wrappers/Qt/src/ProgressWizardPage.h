#ifndef CCPI_PROGRESSWIZARD_PAGE_H
#define CCPI_PROGRESSWIZARD_PAGE_H

#include <QWizardPage>
class ProgressWizardPage :
	public QWizardPage
{
public:
	ProgressWizardPage(void);
	~ProgressWizardPage(void);
	bool isComplete() const;
	void setReconComplete(bool value);
private:
	bool bReconComplete;
};

#endif
