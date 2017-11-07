#include <ui_ReconstructionWizard.h>
#include "ReconstructionImpl.h"
#include "ReconstructionWizardImpl.h"
#include <QFileDialog>
#include <QMessageBox>
#include <QtCore>
#include <QtConcurrent/QtConcurrent>
#include "UserInterfaceQT.h"
extern UserInterfaceQT ui_calls;
ReconstructionWizardImpl::ReconstructionWizardImpl(void)
{
	ui = new Ui_ReconstructionWizard();
	ui->setupUi(this);
	recon = new CCPi::ReconstructionImpl();
	QObject::connect(ui->btnSelectInputFile, SIGNAL(clicked()), this, SLOT(browseInputFile()));
	QObject::connect(ui->btnSelectOutput, SIGNAL(clicked()), this, SLOT(browseOutputFile()));
	QObject::connect(this, SIGNAL(currentIdChanged(int)), this, SLOT(pageChanged(int)));
	QObject::connect(&ui_calls, SIGNAL(send_output(const QString&)), this, SLOT(appendProgressInfo(const QString &)));
	QObject::connect(&ui_calls, SIGNAL(initialise_progress(const int, const QString&)),this,SLOT(initialiseProgressPage(const int,const QString&)));
	QObject::connect(&ui_calls, SIGNAL(update_progress(const int)), ui->pbReconstruction, SLOT(setValue(int)));
	QObject::connect(&watcher, SIGNAL(finished()), this, SLOT(reconstructionComplete()));
	QObject::connect(&saveWatcher, SIGNAL(finished()), this, SLOT(savingComplete()));
	QObject::connect(&ui_calls, SIGNAL(initialise_progress(const int, const QString&)),this,SLOT(initialiseSaveProgressPage(const int,const QString&)));
	QObject::connect(&ui_calls, SIGNAL(update_progress(const int)), ui->pbSavingProgress, SLOT(setValue(int)));
}


ReconstructionWizardImpl::~ReconstructionWizardImpl(void)
{
	delete ui;
	delete recon;
}

bool ReconstructionWizardImpl::validateCurrentPage()
{
	if(this->currentId()==4)
	{
		return SavingOutputPageProcessing();
	}
	return true;
}

void ReconstructionWizardImpl::pageChanged(int id)
{
	switch(id)
	{
	case 1:
		FirstPageProcessing();
		break;
	case 2:
		SecondPageProcessing();
		break;
	case 3:
		ThirdPageProcessing();
		break;
	case 4:
		break;
	}
}

void ReconstructionWizardImpl::FirstPageProcessing()
{
	//Set the device id
	if(ui->btnDevice->isChecked())
		recon->setDeviceId(CCPi::dev_Nikon_XTek);
	//Set the algorithm id
	if(ui->rbAlgorithmCGLS->isChecked())
		recon->setAlgorithmId(CCPi::alg_CGLS);
	if(ui->rbAlgorithmSIRT->isChecked())
		recon->setAlgorithmId(CCPi::alg_SIRT);
	if(ui->rbAlgorithmMLEM->isChecked())
		recon->setAlgorithmId(CCPi::alg_MLEM);
	if(ui->rbAlgorithmTikhonov->isChecked())
		recon->setAlgorithmId(CCPi::alg_CGLS_Tikhonov);
	if(ui->rbAlgorithmTV->isChecked())
		recon->setAlgorithmId(CCPi::alg_CGLS_TVreg);
	//enable hyperthreading
	if(ui->chkBoxHyperThreads->isChecked())
	{
		recon->disableHyperThreads();
	}else{
		recon->enableHyperThreads();
	}
}

void ReconstructionWizardImpl::SecondPageProcessing()
{
	if(ui->leInputFile->text().compare("")!=0)
	{
		recon->setFilename(ui->leInputFile->text().toStdString());
	}
}

void ReconstructionWizardImpl::ThirdPageProcessing()
{
	//Pixels for voxels
	recon->setResolution(ui->spinBoxPixelsPerVoxel->value());
	//number of iterations
	recon->setNumberOfIterations(ui->spinBoxNoOfIterations->value());
	//set beam hardening
	recon->setBeamHardening(ui->chkBoxBeamHardening->isChecked());
	reconstructionProcessId = QtConcurrent::run(recon,&CCPi::ReconstructionImpl::run);
	watcher.setFuture(reconstructionProcessId);
}

bool ReconstructionWizardImpl::SavingOutputPageProcessing()
{
	QString outputfilename = ui->leOutputName->text();
	CCPi::output_format oformat;
	if(ui->rbTiffOutputFormat->isChecked())
		oformat = CCPi::unsigned_short_tiff;
	if(ui->rbRawOutputFormat->isChecked())
		oformat = CCPi::bgs_float_dump;
	//reset the saving is complete
	ui->SavingProgressPage->setSavingComplete(false);
	//Reset the progress bar and text
	ui->lblSavingInfo->setText("");
	ui->pbSavingProgress->setValue(0);

	savingProcessId = QtConcurrent::run(recon,&CCPi::ReconstructionImpl::saveResults, outputfilename.toStdString(), oformat);
	saveWatcher.setFuture(savingProcessId);
	//return recon->saveResults(outputfilename.toStdString(), oformat);
	return true;
}

void ReconstructionWizardImpl::browseInputFile()
{
	if(ui->btnDevice->isChecked()){
		inputFileName = QFileDialog::getOpenFileName(this, tr("Open File"), "", tr("XTek Format(*.xtekct)"));
	}
	ui->leInputFile->setText(inputFileName);
}

void ReconstructionWizardImpl::browseOutputFile()
{
	if(ui->btnGrpOutputFormat->checkedButton()==ui->rbRawOutputFormat)
	{
		outputName = QFileDialog::getSaveFileName(this, tr("Save Raw File"), "", tr("Raw Format(*.raw)"));
	}else if(ui->btnGrpOutputFormat->checkedButton()==ui->rbTiffOutputFormat){
		outputName = QFileDialog::getExistingDirectory(this, tr("Save Tiff File Directory"), "");
	}
	ui->leOutputName->setText(outputName);
}

void ReconstructionWizardImpl::initialiseProgressPage(const int value,const QString& title)
{
	ui->lblProgressTitle->setText(title);
	ui->pbReconstruction->setMaximum(value);
	ui->pbReconstruction->setMinimum(0);
}

void ReconstructionWizardImpl::initialiseSaveProgressPage(const int value,const QString& title)
{
	ui->lblSavingInfo->setText(title);
	ui->pbSavingProgress->setMaximum(value);
	ui->pbSavingProgress->setMinimum(0);
}

void ReconstructionWizardImpl::appendProgressInfo(const QString& info)
{
	ui->tbProgressInfo->append(info);
}

void ReconstructionWizardImpl::reconstructionComplete()
{
	//TODO:Check if the reconstruction thread returned successfully
	ui->ProgressPage->setReconComplete(true);
}

void ReconstructionWizardImpl::savingComplete()
{
	//TODO:Check if the saving thread returned successfully
	ui->SavingProgressPage->setSavingComplete(true);
}