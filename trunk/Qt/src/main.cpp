#include <QApplication>

#include "ReconstructionWizardImpl.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
	ReconstructionWizardImpl *widget = new ReconstructionWizardImpl();
    widget->show();
    return app.exec();
}