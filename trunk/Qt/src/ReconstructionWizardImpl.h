#include <QWizard>
#include <QtCore>
class Ui_ReconstructionWizard;
namespace CCPi{
	class ReconstructionImpl;
};

class ReconstructionWizardImpl :
	public QWizard
{
	Q_OBJECT
public:
	ReconstructionWizardImpl(void);
	~ReconstructionWizardImpl(void);
	bool validateCurrentPage();
public slots:
	void browseInputFile();
	void browseOutputFile();
	void pageChanged(int id);
	void initialiseProgressPage(const int value,const QString& title);
	void initialiseSaveProgressPage(const int value, const QString& title);
	void appendProgressInfo(const QString& info);
	void reconstructionComplete();
	void savingComplete();
private:
	Ui_ReconstructionWizard *ui;
	QString inputFileName;
	QString outputName;
	QFuture<bool> reconstructionProcessId;
	QFuture<bool> savingProcessId;
	QFutureWatcher<bool> watcher;
	QFutureWatcher<bool> saveWatcher;
	CCPi::ReconstructionImpl *recon;
	void FirstPageProcessing();
	void SecondPageProcessing();
	void ThirdPageProcessing();
	bool SavingOutputPageProcessing();
};

